from __future__ import print_function
from Bio import Entrez
import xmltodict
import json
import os
import re
import collections 

os.environ["HOME"] = "/tmp"

def mock_handler(event, context):
    return {
        "statusCode": "200",
        "body": "mockReply"
    }

def build_speechlet_response(title, output, reprompt_text, should_end_session):
    return {
        'outputSpeech': {
            'type': 'PlainText',
            'text': output
        },
        'card': {
            'type': 'Simple',
            'title': "SessionSpeechlet - " + title,
            'content': "SessionSpeechlet - " + output
        },
        'reprompt': {
            'outputSpeech': {
                'type': 'PlainText',
                'text': reprompt_text
            }
        },
        'shouldEndSession': should_end_session
    }


def build_response(session_attributes, speechlet_response):
    return {
        'version': '1.0',
        'sessionAttributes': session_attributes,
        'response': speechlet_response
    }

# --------------- Functions that control the skill's behavior ------------------

def get_welcome_response():
    """ If we wanted to initialize the session to have some attributes we could
    add those here
    """

    session_attributes = {}
    card_title = "Welcome"
    speech_output = "Welcome to the genome wizard by Dark Photon IT consultation for the stars,  " \
                    " I can currently give you a summary of gene function , the location of a gene relative to our chromosomes, and the protein translations of just about any sequenced gene in the human genome , " \
                    " I know about a lot of genes, here is how to learn from me ," \
                    " You can say things like , " \
                    " lookup the gene labelled B R C A 2 , " \
                    " or,  what is the location of T P 5 3 , "  \
                    " or finally , what are the protein translations for G Y F R, " \
                    " Remember , you can use any official gene symbol in the human genome,  which is well over 20 thousand genes available , " \
                    " When reciting the gene symbol , it is best to spell it out letter by letter , " \
                    " Finally , have fun. "
    # If the user either does not reply to the welcome message or says something
    # that is not understood, they will be prompted again with this text.
    reprompt_text = "Please tell me what gene youre looking for by saying for example, " \
                    "tell me about the gene B R C A 2"
    should_end_session = False
    return build_response(session_attributes, build_speechlet_response(
        card_title, speech_output, reprompt_text, should_end_session))


def handle_session_end_request():
    card_title = "Session Ended"
    speech_output = "Thank you for trying the g nome Alexa Skill. " \
                    "Have a nice day! "
    # Setting this to true ends the session and exits the skill.
    should_end_session = True
    return build_response({}, build_speechlet_response(
        card_title, speech_output, None, should_end_session))


def create_gene_label_attribute(gene_label):
    return {"geneLabel": gene_label}


def set_gene_label_in_session(intent, session):
    """ Sets the gene in the session and prepares the speech to reply to the
    user.
    """

    card_title = intent['name']
    session_attributes = {}
    should_end_session = False

    if 'Gene' in intent['slots']:
        geneLabel = intent['slots']['Gene']['value']
        session_attributes = create_gene_label_attribute(geneLabel)
        speech_output = "I now know your selected gene is " + \
                        geneLabel + \
                        ". You can ask me your selected gene by saying, " \
                        "what's my selected gene?"
        reprompt_text = ". You can ask me your selected gene by saying, " \
                        "what's my selected gene??"
    else:
        speech_output = "I'm not sure what your selected gene is. " \
                        "Please try again."
        reprompt_text = "I'm not sure what your selected gene is. " \
                        "You can tell me your selected gene by saying, " \
                        "my selected gene is B R C A 2"
    return build_response(session_attributes, build_speechlet_response(
        card_title, speech_output, reprompt_text, should_end_session))

def prettify_gene_summary(gene_summary):
    result = re.sub('[\[\(].*[\]\)]', '',gene_summary)
    return result 

def prettify_gene_location(gene_label, location):
    if bool(re.match(r"(\d+)(\w{1})(.+)", location)):
        rgx = re.compile('(\d+)(\w{1})(.+)')
        matches = rgx.match(location)
        prettyLoc = "{} is found on the {} arm of chromosome {} at position {}".format(gene_label, matches.group(2), matches.group(1), matches.group(3))
        return prettyLoc
    else:
        prettyLoc = "{} is not on a chromosome but actually at {}".format(gene_label, location)
        return prettyLoc


def get_gene_label_from_session(intent, session):
    session_attributes = {}
    reprompt_text = None

    if session.get('attributes', {}) and "geneLabel" in session.get('attributes', {}):
        selected_gene = session['attributes']['geneLabel']
        speech_output = "Your selected gene is " + selected_gene + \
                        ". Goodbye."
        should_end_session = True
    else:
        speech_output = "I'm not sure what your selected gene is, " \
                        "You can say, my selected gene is B.R.C.A.2."
        should_end_session = False

    # Setting reprompt_text to None signifies that we do not want to reprompt
    # the user. If the user does not respond or says something that is not
    # understood, the session will end.
    return build_response(session_attributes, build_speechlet_response(
        intent['name'], speech_output, reprompt_text, should_end_session))


def gene_encodes_proteins(library):
    if "Entrezgene_type" in library["Entrezgene-Set"]["Entrezgene"]:
        geneType = library["Entrezgene-Set"]["Entrezgene"]["Entrezgene_type"].items()[0][1]
        proteinCount = library["Entrezgene-Set"]["Entrezgene"]["Entrezgene_type"].items()[1][1]
        if geneType == 'protein-coding' and proteinCount > 0:
            return True
        else:
            return False              

def get_proteins_encoded_by_gene(library):
    proteins = list()
    if gene_encodes_proteins(library):
        geneType = library["Entrezgene-Set"]["Entrezgene"]["Entrezgene_type"].items()[0][1]
        proteinCount = library["Entrezgene-Set"]["Entrezgene"]["Entrezgene_type"].items()[1][1]
        if "Entrezgene_prot" in library["Entrezgene-Set"]["Entrezgene"]:
            if "Prot-ref" in  library["Entrezgene-Set"]["Entrezgene"]["Entrezgene_prot"]:
                if "Prot-ref_name" in library["Entrezgene-Set"]["Entrezgene"]["Entrezgene_prot"]["Prot-ref"]:
                    for k in library["Entrezgene-Set"]["Entrezgene"]["Entrezgene_prot"]["Prot-ref"]["Prot-ref_name"].keys():
                        proteins.extend(library["Entrezgene-Set"]["Entrezgene"]["Entrezgene_prot"]["Prot-ref"]["Prot-ref_name"][k])
    return geneType, proteinCount, proteins

def prettify_protein_list(gene_label,proteins,protein_count):
    speech_output = str()
    protein_count = int(protein_count)
    speech_output += "{} translates for {} proteins  , here is a list of those proteins  , ".format(gene_label, protein_count)
    cnt = 0
    for protein in proteins:
        cnt += 1
        if (cnt-1) < protein_count:
            speech_output += " {} ,".format(protein)
        else:
            speech_output += " and finally {} .".format(protein)
    return speech_output              

def query_entrez_proteins(gene_label):
    session_attributes = {"geneLabel": gene_label}
    card_title = "ProteinsTranslatedFor"

    print("Setting OS Environ Values")
    Entrez.email = os.environ.get("ENTREZ_EMAIL")
    Entrez.api_key = os.environ.get("ENTREZ_API_KEY")
    print("Gathering records for Homo Sapiens {} gene".format(gene_label))
    handle = Entrez.esearch(db="gene", term="(\"Homo Sapiens\"[orgn]) AND \"%s\"[Gene]"%(gene_label)) #pass in query terms to search.
    record = Entrez.read(handle)
    if int(record['Count']) == 0:
        print("Query Returned No Results, Try spelling out the official HGNC gene symbol")
        speech_output = "Query Returned No Results, Try spelling out the official HGNC gene symbol as your search string. Visit www.genenames.org and search by keyword or symbol for references."
        reprompt_text = None
        should_end_session = False
        print(build_response(session_attributes, build_speechlet_response(card_title, speech_output, reprompt_text, should_end_session)))
        return build_response(session_attributes, build_speechlet_response(card_title, speech_output, reprompt_text, should_end_session))
    elif int(record['Count']) == 1:
        print("Query was successful, here are the results.")
        handle = Entrez.efetch(db="gene", id=record['IdList'][0], retmode="xml")
        library = xmltodict.parse(handle)
        handle.close()
        #tmploc = str(library["Entrezgene-Set"]["Entrezgene"]["Entrezgene_location"]["Maps"]["Maps_display-str"])
        if gene_encodes_proteins(library):    
            geneType, proteinCount, proteinList = get_proteins_encoded_by_gene(library)
            speech_output =  prettify_protein_list(gene_label , proteinList, proteinCount)
        else:
            speech_output = "This gene does not translate to any proteins."
        #speech_output = str(library["Entrezgene-Set"]["Entrezgene"]["Entrezgene_location"]["Maps"]["Maps_display-str"])
        reprompt_text = None
        should_end_session = False
        print(build_response(session_attributes, build_speechlet_response(card_title, speech_output, reprompt_text, should_end_session)))
        return build_response(session_attributes, build_speechlet_response(card_title, speech_output, reprompt_text, should_end_session))
    elif int(record['Count']) > 1:
        print("Query Returned Multiple Results, Try spelling out the gene symbol as officially defined by HGNC, or be more specific about which gene you would like to know about.")
        speech_output = "Query Returned Multiple Results, Try spelling out the gene symbol as officially defined by HGNC, or be more specific about which gene you would like to know about."
        reprompt_text = None
        should_end_session = False
        print(build_response(session_attributes, build_speechlet_response(card_title, speech_output, reprompt_text, should_end_session)))
        return build_response(session_attributes, build_speechlet_response(card_title, speech_output, reprompt_text, should_end_session))

def query_entrez_location(gene_label):
    session_attributes = {"geneLabel": gene_label}
    card_title = "LocateGene"

    print("Setting OS Environ Values")
    Entrez.email = os.environ.get("ENTREZ_EMAIL")
    Entrez.api_key = os.environ.get("ENTREZ_API_KEY")
    print("Gathering records for Homo Sapiens BRCA2 gene")
    handle = Entrez.esearch(db="gene", term="(\"Homo Sapiens\"[orgn]) AND \"%s\"[Gene]"%(gene_label)) #pass in query terms to search.
    record = Entrez.read(handle)
    if int(record['Count']) == 0:
        print("Query Returned No Results, Try spelling out the official HGNC gene symbol")
        speech_output = "Query Returned No Results, Try spelling out the official HGNC gene symbol as your search string. Visit www.genenames.org and search by keyword or symbol for references."
        reprompt_text = None
        should_end_session = False
        print(build_response(session_attributes, build_speechlet_response(card_title, speech_output, reprompt_text, should_end_session)))
        return build_response(session_attributes, build_speechlet_response(card_title, speech_output, reprompt_text, should_end_session))
    elif int(record['Count']) == 1:
        print("Query was successful, here are the results.")
        handle = Entrez.efetch(db="gene", id=record['IdList'][0], retmode="xml")
        library = xmltodict.parse(handle)
        handle.close()
        tmploc = str(library["Entrezgene-Set"]["Entrezgene"]["Entrezgene_location"]["Maps"]["Maps_display-str"])
        speech_output = prettify_gene_location(gene_label, tmploc)
        reprompt_text = None
        should_end_session = False
        print(build_response(session_attributes, build_speechlet_response(card_title, speech_output, reprompt_text, should_end_session)))
        return build_response(session_attributes, build_speechlet_response(card_title, speech_output, reprompt_text, should_end_session))
    elif int(record['Count']) > 1:
        print("Query Returned Multiple Results, Try spelling out the gene symbol as officially defined by HGNC, or be more specific about which gene you would like to know about.")
        speech_output = "Query Returned Multiple Results, Try spelling out the gene symbol as officially defined by HGNC, or be more specific about which gene you would like to know about."
        reprompt_text = None
        should_end_session = False
        print(build_response(session_attributes, build_speechlet_response(card_title, speech_output, reprompt_text, should_end_session)))
        return build_response(session_attributes, build_speechlet_response(card_title, speech_output, reprompt_text, should_end_session))

def query_entrez_summary(gene_label):
    session_attributes = {"geneLabel": gene_label}
    card_title = "DescribeGene"
    
    print("Setting OS Environ Values")
    Entrez.email = os.environ.get("ENTREZ_EMAIL")
    Entrez.api_key = os.environ.get("ENTREZ_API_KEY")
    print("Gathering records for the HomoSapiens {} gene".format(gene_label))
    handle = Entrez.esearch(db="gene", term="(\"Homo Sapiens\"[orgn]) AND \"%s\"[Gene]"%(gene_label)) #pass in query terms to search.
    record = Entrez.read(handle)


    if int(record['Count']) == 0:
        print("Query Returned No Results, Try spelling out the official HGNC gene symbol")
        speech_output = "Query Returned No Results, Try spelling out the official HGNC gene symbol as your search string. Visit www.genenames.org and search by keyword or symbol for references."
        reprompt_text = None
        should_end_session = False
        print(build_response(session_attributes, build_speechlet_response(card_title, speech_output, reprompt_text, should_end_session)))
        return build_response(session_attributes, build_speechlet_response(card_title, speech_output, reprompt_text, should_end_session))
    elif int(record['Count']) == 1:
        print("Query was successful, here are the results.")
        handle = Entrez.efetch(db="gene", id=record['IdList'][0], retmode="xml")
        library = xmltodict.parse(handle)
        handle.close()
        speech_output = prettify_gene_summary(library["Entrezgene-Set"]["Entrezgene"]["Entrezgene_summary"])
        reprompt_text = None
        should_end_session = False
        print(build_response(session_attributes, build_speechlet_response(card_title, speech_output, reprompt_text, should_end_session)))
        return build_response(session_attributes, build_speechlet_response(card_title, speech_output, reprompt_text, should_end_session))
    elif int(record['Count']) > 1:
        print("Query Returned Multiple Results, Try spelling out the gene symbol as officially defined by HGNC, or be more specific about which gene you would like to know about.")
        speech_output = "Query Returned Multiple Results, Try spelling out the gene symbol as officially defined by HGNC, or be more specific about which gene you would like to know about."
        reprompt_text = None
        should_end_session = False
        print(build_response(session_attributes, build_speechlet_response(card_title, speech_output, reprompt_text, should_end_session)))
        return build_response(session_attributes, build_speechlet_response(card_title, speech_output, reprompt_text, should_end_session))

# --------------- Events ------------------

def on_session_started(session_started_request, session):
    """ Called when the session starts """

    print("on_session_started requestId=" + session_started_request['requestId']
          + ", sessionId=" + session['sessionId'])


def on_launch(launch_request, session):
    """ Called when the user launches the skill without specifying what they
    want
    """

    print("on_launch requestId=" + launch_request['requestId'] +
          ", sessionId=" + session['sessionId'])
    # Dispatch to your skill's launch
    return get_welcome_response()


def on_intent(intent_request, session):
    """ Called when the user specifies an intent for this skill """

    print("on_intent requestId=" + intent_request['requestId'] +
          ", sessionId=" + session['sessionId'], ", intentName=" + intent_request['intent']['name'])

    intent = intent_request['intent']
    intent_name = intent_request['intent']['name']
    #print("Gene Selected: %s"%(intent['slots']['Gene']['value']))

    # Dispatch to your skill's intent handlers
    if intent_name == "DescribeGene":
        return query_entrez_summary(str(intent['slots']['Gene']['value']))
    elif intent_name == "LocateGene":
        return query_entrez_location(str(intent['slots']['Gene']['value']))
    elif intent_name == "ProteinsTranslatedFor":
        return query_entrez_proteins(str(intent['slots']['Gene']['value']))
    elif intent_name == "AMAZON.HelpIntent":
        return get_welcome_response()
    elif intent_name == "AMAZON.CancelIntent" or intent_name == "AMAZON.StopIntent":
        return handle_session_end_request()
    else:
        raise ValueError("Invalid intent")


def on_session_ended(session_ended_request, session):
    """ Called when the user ends the session.

    Is not called when the skill returns should_end_session=true
    """
    print("on_session_ended requestId=" + session_ended_request['requestId'] +
          ", sessionId=" + session['sessionId'])
    # add cleanup logic here





# --------------- Main handler ------------------

def lambda_handler(event, context):
    #print(event.keys())

    if event['session']['new']:
        on_session_started({'requestId': event['request']['requestId']},event['session'])
    if event['request']['type'] == "LaunchRequest":
        return on_launch(event['request'], event['session'])
    elif event['request']['type'] == "IntentRequest":
        return on_intent(event['request'], event['session'])
    elif event['request']['type'] == "SessionEndedRequest":
        return on_session_ended(event['request'], event['session'])

