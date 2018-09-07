import nltk
import os
import csv
import numpy as np
import math

# TODO update these paths
# export PYTHONPATH=~/modules
# export CLASSPATH=~/Teaching/C4M/Code/stanford-parser-full-2015-12-09/

from nltk.stem.porter import *
from nltk.tag.stanford import StanfordPOSTagger
from nltk.internals import find_jars_within_path
from nltk.parse.stanford import StanfordParser
from sklearn import svm


def load(file):
    """ (file open for reading) -> list of str

    Return the list of filenames from file, which contains
    one filename per line.
    """

    file_list = []
    for address in file:
        file_list.append(address.strip('\n'))
    return file_list


def preprocess(flist, folder_path):
    """ (file open for reading, str) -> Nonetype

    flist contains one filename per line and folder_path represents a 
    directory. Do preprocessing on each file from flist in folder_path.
    """

    error_log = []
    for i in range(len(flist)):

        path = flist[i]

        stemmer = PorterStemmer()
        parser = StanfordParser(
            model_path='edu/stanford/nlp/models/lexparser/englishPCFG.ser.gz',
             verbose=True)
        stanford_dir = parser._classpath[0].rpartition('/')[0]
        parser._classpath = tuple(find_jars_within_path(stanford_dir))

        with open(path, 'r') as rf:
            try:
                sent = [line.strip('\n ') for line in rf]
            except UnicodeDecodeError as e:
                error_log.append('Unicode Decode Error:\t' + path + '\n')
                pass
            else:
                if not sent:
                    error_log.append('Empty File Error:\t' + path + '\n')
                    pass
                else:
                    # Stemming with Porter Stemmer
                    pars_stem = stemmer.stem(' '.join(sent))
                    stemmed = '\n'.join(sent)

                    wf = open(folder_path 
                        + path.split('.')[0].split('/')[-1] + '.stem', 'w')
                    wf.write(stemmed)
                    wf.close()

                    # POS Tagging after tokenizing and stemming
                    pos = nltk.pos_tag(pars_stem.split())
                    wf = open(folder_path 
                        + path.split('.')[0].split('/')[-1] + '.pos', 'w')
                    wf.write(str(pos))
                    wf.close()

                    # CFG parser
                    try:
                        parsed = parser.raw_parse(pars_stem)
                    except (TypeError, IndexError, NameError) as e:
                        error_log.append('Unparsable Error:/t' + path + '/n')
                        pass
                    wf = open(folder_path 
                        + path.split('.')[0].split('/')[-1] + '.pars', 'w')
                    s_pars = " ".join(str(x) for x in list(parsed))
                    s_pars = s_pars.replace("Tree", "")
                    s_pars = s_pars.replace("[","")
                    s_pars = s_pars.replace("]","")
                    s_pars = s_pars.replace("\'","")
                    wf.write(s_pars)
                    wf.close()

    # Print files paths with Errors
    if error_log:
        wf = open(folder_path + 'error_log', 'wb')
        for line in error_log:
            wf.write(line)
        wf.close()


def read_norms(file):
    """ (file open for reading) -> dict of {str: list of str}

    Read the contents of file and return a dictionary where each key is a
    word and each value is a list of strings representing AoA, IMG, and FAM.

    """

    norms = {}
    with open(file, 'r') as csvfile:
        myreader = csv.reader(csvfile, delimiter=',')
        next(myreader) # skip header
        for row in myreader:
            norms[ row[1]] = [row[3], row[4], row[5] ] # Get AoA, IMG, FAM

    return norms


# #### TODO (BELOW) ####
# Add helper functions to computer feature values.

# Feature 1
def count_words(stem_file):
    """ (stem file open for reading) -> int

    Return the number of words in stem_file.
    """

    count = 0
    for line in stem_file:
        count += 1
    return count

# #### TODONE (ABOVE) #### 

# Feature 2
def avg_length(f):
    """ (stem file open for reading) -> float

    Return the average number of characters in the utterance in f.
    """

    total = 0
    count = 0
    for line in f:
        total += len(line.strip()) 
        count += 1
    return total / count

# Feature 3: 
def honore(f):
    """ (stem file open for reading) -> float

    Return the Honore statistic.
    """
    
    word_to_count = {}
    for word in f:
        word_to_count[word.strip()] = word_to_count.get(word.strip(), 0) + 1
    num_tokens = sum(list(word_to_count.values()))
    num_types = len(word_to_count)
    num_types_once = (list(word_to_count.values())).count(1)
    
    # avoiding division by 0
    honore = 0
    if num_types_once != num_types:
        honore = (100 * math.log(num_tokens)) / (1 - (num_types_once / num_types)) 
    return honore

def parse_tree_depth(f):
    """ (pars file open for reading) -> float

    Return the parse tree depth.
    """

    utterance = f.read()
    max_depth = 0
    cur_depth = 0
    for ch in utterance:
        if ch == ')':
            cur_depth += 1
        else:
            max_depth = max(max_depth, cur_depth)
            cur_depth = 0
    return max_depth

# Feature 7
def aoa(f, norms):
    """ (stem file open for reading) -> float

    Return the AoA.
    """
    total = 0
    count = 0
    for line in f:
        word = line.strip()
        if word in norms:
            result = norms[word][0]
            total += int(norms[word][0])
            count += 1

    # avoiding division by 0
    avg = 0
    if count != 0:
        avg = total / count
    return avg

def feature_eight(f):
    """ (pos file open for reading) -> float

    Return the result of calculating feature eight.
    """
    parse = f.read()
    NN = parse.count('NN')
    NNS = parse.count('NNS')
    NNP = parse.count('NNP')
    NNPS = parse.count('NNPS')
    PRP = parse.count('PRP')
    PRP_2 = parse.count('PRP$')

    # avoiding division by 0
    result = 0
    if (PRP + PRP_2):
        result = (NN + NNS + NNP + NNPS) / (PRP + PRP_2)
    return result

def extract_features(flist, path):
    """ (file open for reading, str) -> array

    Return an N x 9 array with nine features for each file in flist
    (one filename per line) in directory path.
    """

    norms = read_norms('Norms.csv')
    features = np.array([]) 

    for filename in flist:

        # #### TODO (BELOW) ####
        # 0. Count number of words in utterance
        f = open(path + filename + '.stem', 'r')
        f0 = count_words(f)
        f.close()

        # 1. Count average number of characters in utterance
        f = open(path + filename + '.stem', 'r')
        f1 = avg_length(f)
        f.close()

        # 2. Compute Honore's statistic on utterance
        f = open(path + filename + '.stem', 'r')
        f2 = honore(f)
        f.close()
  
        # 3. Compute the parse tree depth
        f = open(path + filename + '.pars', 'r')
        f3 = parse_tree_depth(f)
        f.close()

        # 4. Count the number of 'CC' instances in parse
        f = open(path + filename + '.pos', 'r')
        parse = f.read()
        count = parse.count('CC')
        f4 = count
        f.close()

        # 5. Count the number of 'VBG' instances in parse
        f = open(path + filename + '.pos', 'r')
        parse = f.read()
        count = parse.count('VBG')
        f5 = count
        f.close()

        # 6. Count the number of 'VBZ' and 'VBP' instances in parse
        f = open(path + filename + '.pos', 'r')
        parse = f.read()
        count = parse.count('VBZ') + parse.count('VBP')
        f6 = count
        f.close()

        # 7. Count the average Age of Acquisiton (AoA) of words
        f = open(path + filename + '.stem', 'r')
        f7 = aoa(f, norms)
        f.close()

        # 8. Compute (NN + NNS + NNP + NNPS) / (PRP + PRP$)
        f = open(path + filename + '.pos', 'r')
        f8 = feature_eight(f)
        f.close()

        vector = np.array([f0, f1, f2, f3, f4, f5, f6, f7, f8])
        features = np.vstack([features, vector]) if features.size else vector 

        # #### TODONE (ABOVE) #### 
    
    return features


def classify(DD, CD):
    """ (array, array) -> Nonetype

    Report the outcome of classifying feature vectors.  DD and CD are arrays for 
    participants with and without dementia, respectively.
    """

    # do K-fold cross validation
    all_data  = np.concatenate((DD,CD))
    all_class = np.vstack((np.ones((DD.shape[0], 1)), np.zeros((CD.shape[0], 1))))
    N = all_data.shape[0]
    K = 5
    accuracies = np.zeros((K, 1))
    sensitivities = np.zeros((K, 1))
    specificities = np.zeros((K, 1))
    randIndices = np.random.permutation(range(0, N))

    for fold in range(0, K) :
        i_test = randIndices[fold * (N // K) : (fold + 1) * (N // K)]
        i_train = [val for val in randIndices if val not in i_test]
        c_train = all_class[i_train]
        c_test = all_class[i_test]   

        # train the model using features all_data[i_train,:] and classes c_train 
        # TODO
        my_model = svm.SVC(kernel="linear")
        c_train = np.ravel( c_train )
        my_model.fit( all_data[i_train, :] , c_train )
        
        my_predictions = my_model.predict(all_data[i_test, :])

        # # compute the % errors using all_data[i_test,:], c_test, and the model you trained
        # # TODO
        TP = len( [i for i in range(0,c_test.size) if my_predictions[i]==1 and c_test[i]==1] )
        TN = len( [i for i in range(0,c_test.size) if my_predictions[i]==0 and c_test[i]==0] )
        FP = len( [i for i in range(0,c_test.size) if my_predictions[i]==1 and c_test[i]==0] )
        FN = len( [i for i in range(0,c_test.size) if my_predictions[i]==0 and c_test[i]==1] )

        accuracies[fold] = (TP + TN) / (TP + TN + FP + FN)
        sensitivities[fold] = TP / (TP + FN)
        specificities[fold] = TN / (TN + FP)


    # report result mean and variance, over all folds, of each of accuracy, specificity, and sensitivity
    # TODO
    print('mean accurary', np.mean(accuracies))
    print('mean sensitivity', np.mean(sensitivities))
    print('mean specificity', np.mean(specificities))
    print('variance accurary', np.var(accuracies))
    print('variance sensitivity', np.var(sensitivities))
    print('variance specificity', np.var(specificities))
    

if __name__ == "__main__":

    # read list of txt2 files
    df = open("./results/Dementia.list", "r")
    cf = open("./results/Controls.list", "r")

    # load list of txt2 files to be parsed
    dlist = load(df)
    clist = load(cf)

    # close the files
    df.close()
    cf.close()

    # do preprocessing. (We've already taken care of this.)
    #preprocess(dlist, "./results/Dementia/")
    #preprocess(clist, "./results/Controls/")
    
    # extract relevant features
    DD = extract_features(dlist, "./Data/Dementia/")
    CD = extract_features(clist, "./Data/Controls/")

    # do the K-fold crossvalidation classification
    classify(DD, CD )

