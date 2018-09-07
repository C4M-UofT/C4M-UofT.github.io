import os
import math

HOMEFOLDER = "/media/guerzhoy/Windows/cs_for_docs/workshops/diagnosis_tfidf/wikipages/"    

def replace_punctuation_w_space(text):
    punctuation = ['"', "'", "<", ">", ",", "!", "?", ".", "\n"]
    for p in punctuation:
        text = text.replace(p, " ")
    
    return text


def clean_up_file(text):
    text = text.lower()
    text = replace_punctuation_w_space(text)
    return text
    

def get_all_texts(folder):
    filenames = os.listdir(folder)
    
    
    all_texts = {}
    
    for filename in filenames:
        text = open(folder + filename).read()
        text = clean_up_file(text)
        disease = filename.split(".")[0]
        all_texts[disease] = text
    
    return all_texts
    
def log_inv_doc_freq(all_texts, term):
    count = 0
    for disease, text in all_texts.items():
        if text.count(term) > 0:
            count += 1
    
    if count == 0:
        return -1
    else:
        return math.log(len(all_texts)/count)


def add_weighted_term_freq(all_texts, term, scores):
    idf = log_inv_doc_freq(all_texts, term)
    if idf < 0:
        return
    for disease, text in all_texts.items():
        if term in text:
            scores[disease] += idf
            
    
def compute_scores(all_texts, terms):
    scores = {}
    for disease in all_texts:
        scores[disease] = 0
    
    for term in terms:
        add_weighted_term_freq(all_texts, " " + term + " ", scores)
        
    return scores
    
def best_match(scores):
    k = list(scores.keys())
    best_score = scores[k[0]]
    disease = k[0]
    
    for d, score in scores.items():
        if score > best_score:
            best_score = score
            disease = d
            
    return disease

def query(folder, terms):
    all_texts = get_all_texts(folder)
    scores = compute_scores(all_texts, terms)    
    return best_match(scores)



while True:
    q = clean_up_file(input("query: ")).split()
    print("Best match:", query(HOMEFOLDER, q))

