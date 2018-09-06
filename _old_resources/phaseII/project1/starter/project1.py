import os
import math

# Redefine this constant to be the path to where you stored the unzipped 
# wikipages. The Windows OS uses \ as the directory separator and inside
# this constant string you need two slashes for each separator.
# For example "C:\\user\\folder\\subfolder\\"
HOMEFOLDER = "wikipages"    
PUNCTUATION = """.,<>;'":{}[]|!@#$%^&*()"""

def clean_up(text):
    """ (str) -> str

    Return a new version of text in which all the letters have been 
    converted to lowercase, and all punctuation is replaced with spaces.

    >>> clean_up('Influenza, commonly known as "the flu", is ...')
    'influenza  commonly known as  the flu   is    '
    """

def get_documents(datapath):
    """ (str) -> dict of {str: str}

    Return a dictionary where the keys are document names (without .html)
    and the values are the contents of the file corresponding .html file
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
            doc_to_text[disease] = text
    
    return doc_to_text
    

def keyword_found(keyword, doc_name, doc_to_text):
    """ (str, str, dict of {str:str}) -> bool
    
    Return True iff keyword is found in this doc_name inside doc_to_text
    as a full token separated by whitespace.
    
    """

    
def idf(keyword, doc_to_text):
    """ (str, dict of {str: str}) -> float

    Return the idf for this keyword in documents doc_to_text. 

    """


def build_empty_scores_dict(doc_to_text):
    """ (dict of {str: str}) -> dict of {str: number}

    Build and return an empty dictionary where the keys are the same as the keys in doc_to_text
    and the values are all 0.
    """
        

def update_scores(doc_to_score, keyword, doc_to_text):
    """ (dict of {str: number}, str, dict of {str: str}) -> None

    Update current_scores by adding to the value of each entry to TF-IDF individual score
    for keyword based on the documents in all_texts.
    """
    
