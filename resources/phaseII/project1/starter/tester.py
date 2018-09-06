# Here are a few very simple tests to run on your functions to ensure that they
# have the signatures that we specified. (A signature of a function is the
# type and order of the parameters and also the type of the returned value.)
#
# These basic tests do NOT thoroughly test your code. Do that yourself.
# They provide a confirmation that your understanding of the types of each
# parameter was correct. If one of your functions does not pass one of these 
# tests or does not even complete, then something is wrong and you should 
# come to office hours before submitting your solution.
#
# Only the first two tests specify the output. For the next five, you'll need
# to read the input and work out the output that you should be expecting.
#
################################################################################
#Testing keyword_found
#
#sample dictionary
doc_to_text = {"doc1": "mykeywords mykeywords", "doc2": "mykeyword mykeyword"}

actual   = keyword_found("mykeyword", "doc1", doc_to_text)
expected = False 
print("TEST 1 output:", actual, "\tTEST 1 expected:", expected, "\tPASSED?", actual == expected)


actual   = keyword_found("mykeyword", "doc2", doc_to_text)
expected = True 
print("TEST 2 output:", actual, "\tTEST 2 expected:", expected, "\tPASSED?", actual == expected)



################################################################################
#Test idf

#sample dictionary
doc_to_text = {"doc1": "a b c", "doc2": "c d f", "doc3": "g e h"}

actual = idf("z", doc_to_text)
print("TEST 3 output:", actual)

actual = idf("a", doc_to_text)
print("TEST 4 output:", actual)

actual = idf("c", doc_to_text)
print("TEST 5 output:", actual)

################################################################################
#Test build_empty_scores_dict
doc_to_text = {"doc1": "a b c", "doc2": "c d f", "doc3": "g e h"}
actual = build_empty_scores_dict(doc_to_text)
print("TEST 6 output:", actual)

################################################################################
#Test update_scores

doc_to_text = {"doc1": "a b c", "doc2": "c d f", "doc3": "g e h"}
doc_to_score = {"doc1": 0, "doc2": 0, "doc3": 0}
update_scores(doc_to_score, "c", doc_to_text)
print("TEST 7 output:", doc_to_score)

