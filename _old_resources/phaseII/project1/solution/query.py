import os
import math

# Redefine this constant to be the path to where you stored the unzipped 
# wikipages. The Windows OS uses \ as the directory separator and inside
# this constant string you need two slashes for each separator.
# For example "C:\\user\\folder\\subfolder\\"
#HOMEFOLDER = "/home/guerzhoy/Desktop/diagnosis/wikipages/"    
#HOMEFOLDER = "/Users/mcraig/admin/Activity2016/cs_for_docs/workshops/diagnosis_tfidf/wikipages/"
HOMEFOLDER = '/Users/campbell/courses/cs_for_docs/cohort2_sep2016/phase2/project1/wikipages/'
PUNCTUATION = """.,<>;'":{}[]|!@#$%^&*()"""

def clean_up(text):
    """ (str) -> str

    Return a version of the string text where all the letters have been 
    converted to lower case, and all punctuation was replaced with whitespaces

    >>> clean_up('Influenza, commonly known as "the flu", is ...')
    'influenza  commonly known as  the flu   is    '
    """
    result = ""
    for ch in text:
        if ch in PUNCTUATION:
            result += " "
        else:
            result += ch
    return result
    


def get_documents(datapath):
    """ (str) -> dict of {str: str}

    Return a dictionary where the keys are disease names
    and the values are the contents of the file key.html
    from the directory datapath.
    """
    
    # get a list of all the filenames in the directory
    filenames = os.listdir(datapath)
    
    # dictionary of all texts, keys are disease names
    doc_to_text= {}
    
    for filename in filenames:

        # only consider filenames that end in ".html"
        if len(filename) > 5  and  filename[-5:] == ".html":

            # read the entire file's contents as a string
            text = open(datapath + filename).read()
            # since all the filenames end in .html, just drop that part
            disease = filename[:-5]
            # insert it into the dictionary
            doc_to_text[disease] = clean_up(text)
    
    return doc_to_text
    

def keyword_found(keyword, doc, doc_to_text):
    """ (str, str, dict of {str:str}) -> bool
    
    Return True iff keyword is found in doc_to_text[doc] as a full
    token separated by whitespace.
    
    """
    text = doc_to_text[doc]
    tokens = text.split()
    return keyword in tokens

    
def idf(keyword, doc_to_text):
    """ (str, dict of {str: str}) -> float

    Return the idf for this keyword in this collection of documents. If the keyword
    is in 0 documents, return -1.
    
    """
    docs_containing_keyword = 0
    for disease in doc_to_text:
        if keyword_found(keyword, disease, doc_to_text):
            docs_containing_keyword += 1
            print(disease)
    #print("docs_containing", keyword, docs_containing_keyword)
    if docs_containing_keyword == 0:
        return -1
    print('num documents', len(doc_to_text))
    return math.log(len(doc_to_text)/docs_containing_keyword)


def build_empty_scores_dict(doc_to_text):
    """ (dict of {str:str}) -> dict of {str:number}
    Build and return an empty dictionary where the keys are the same as the keys in doc_to_text
    and the values are all 0.
    """
    result = {}
    for key in doc_to_text:
        result[key] = 0.0
    return result
        
def update_scores(doc_to_score, keyword, doc_to_text):
    """ (dict of {str: number}, str, dict of {str: str}) -> None

    Update current_scores by adding to the value of each entry to TF-IDF individual score
    for keyword based on the documents in doc_to_text.

    """
    # since this value is the same for all diseases, just compute it once
    shared_idf = idf(keyword, doc_to_text)

    # now update the dictionary for each entry
    for disease in doc_to_score:
        if keyword_found(keyword, disease, doc_to_text):
            doc_to_score[disease] += shared_idf
        # don't need an else clause since if keyword is not found score stays the same
 
    # no need to return any value since dictionary passed as parameter was changed


# before you can call this function, you need to write functions for 
#   clean_up
#   query  
# and define global constant HOMEFOLDER
def interface():
    while True:
        q = clean_up(input("query: ")).split()
        print("Best match:", query(HOMEFOLDER, q))


def testing_scores(scores):
    print("Mumps", scores["Mumps"], "\tTyphus", scores["Typhus"], "\tMeningitis", scores["Meningitis"])


texts = get_documents(HOMEFOLDER)

#print("fever Typhus ",tf("fever", "Typhus", texts))
#print(tf("fever", "Mumps", texts))
#print("deafness Mumps", tf("deafness", "Mumps", texts))
#print("deafness Meningitis ", tf("deafness", "Meningitis", texts))
#print("deafness Typhus ", tf("deafness", "Typhus", texts))
#print("bacteria Mumps", tf("bacteria", "Mumps", texts))
#print("bacteria Meningitis ", tf("bacteria", "Meningitis", texts))
#print("bacteria Typhus ", tf("bacteria", "Typhus", texts))
#
#print("idf for fever", idf("fever", texts))
#print("idf for bacteria", idf("bacteria", texts))
#print("idf for deafness", idf("bacteria", texts))
#print("idf for the", idf("the", texts))
#print("idf for parotid", idf("parotid", texts))
#
#scores = build_empty_scores_dict(texts)
#testing_scores(scores)
#update_scores(scores, "parotid", texts)
#testing_scores(scores)
#update_scores(scores, "the", texts)
#testing_scores(scores)
#update_scores(scores, "bacteria", texts)
#testing_scores(scores)
#update_scores(scores, "deafness", texts)
#testing_scores(scores)

#print(idf("notfound", texts))

query = input("enter a list of keywords separated by spaces or 'quit' to stop: ")
while query != "quit":
    keywords = clean_up(query).split()
    scores = build_empty_scores_dict(texts)
    for keyword in keywords:
        update_scores(scores, keyword, texts)

    # now find best match
    best_score = 0
    current_best_disease = "no match"
    for disease, score in scores.items():
        if score > best_score:
            current_best_disease = disease
            best_score = score
    print(best_score)
    print("The best match was ", current_best_disease)
    query = input("enter a list of keywords separated by spaces or 'quit' to stop: ")
    
    
