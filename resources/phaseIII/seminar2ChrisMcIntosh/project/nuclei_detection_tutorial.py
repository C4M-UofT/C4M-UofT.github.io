import numpy as npimport pylabimport sysfrom skimage import iofrom skimage import colorfrom scipy import ndimagefrom skimage.feature import match_templatefrom skimage.feature import peak_local_maxfrom matplotlib import pyplotimport project_helpers as ph# Flag that indicates whether to show images.VIEW = False# Set the path to the data files.data_path = "crchistophenotypes_2016_04_28/CRCHistoPhenotypes_2016_04_28/"#################Nucleus detection - Template matching Tutorial##################Click the pan-axis button, and pan while holding left mouse; zoom while holding right#Find the boundary box of a template cell, or use the one belowpyplot.imshow(img)pyplot.autoscale(False)pyplot.plot(detections[:, 0], detections[:, 1], 'kx')if(VIEW):    pyplot.show()#Create the template nucleusT = img_gray[255:265,178:192]pylab.imshow(T)if(VIEW):    pyplot.show()#Use template matching, which is very similar to convolution#to find out where the template appears in the image.match = match_template(img_gray, T, pad_input=True)pylab.imshow(match)if(VIEW):    pyplot.show()#Use non-maximal supression to find peaks in the output that#indicate our template had a strong match at that positioncoordinates = peak_local_max(match, min_distance=20)coordinates = np.fliplr(coordinates) #Flip the results to plot(x,y) from matrix (row,column)###################TODO 2.1 finish the score_detector function to score the results###################precision, recall, F1Measure, FNList = ph.score_detector(detections, coordinates)#Display the results- There are a lot of false negatives!ph.display_detection_results(img, detections, coordinates, FNList, VIEW)#################Nucleus detection - Features and Machine learning##################Load the image and the nuclei dataimg, detections, labels = ph.load_data_set(data_path,"img7")#Get the same template used earlierimg_gray = ph.my_rgb2gray(img)T = img_gray[255:265, 178:192]#Part 3.1 - Following the tutorial, fill in calculate_features in project_helpers.pyfeaturesND = ph.calculate_features(img,T)#Get some training samplesN = 15000features_train, labels_train = ph.get_training_samples(featuresND, detections, img_gray, N)        classif_oneTemplate, scaler_oneTemplate = ph.train_model(features_train, labels_train)#Predict the outputprobability_prediction, predictions = ph.predict_centres(img_gray, featuresND, scaler_oneTemplate, classif_oneTemplate)#Display the resulting probability distribution showing where it thinks centres arepyplot.imshow(probability_prediction)if(VIEW):    pyplot.show()#Score the new results and compare to the imageProcessing_Tutorial.py - A huge jump in true positives!precision, recall, F1Measure, FNList = ph.score_detector(detections, predictions)ph.display_detection_results(img, detections, predictions, FNList,VIEW)#Part 3.2 - What about the nuclei that look different from our template?#Get an entire set of templates! One for each type of celltemplates = ph.get_templates(img_gray, detections, 25)#Visualize the set of templatesph.display_templates(templates, VIEW)#Train a modelfeaturesND = ph.calculate_features(img, templates)features_train,labels_train = ph.get_training_samples(featuresND, detections, img_gray, N)classif_multipleTemplates, scaler_multipleTemplates = ph.train_model(features_train, labels_train)        #Predict the outputprobability_prediction, predictionsNew = ph.predict_centres(img_gray, featuresND, scaler_multipleTemplates, classif_multipleTemplates)#Score the old resultsprecision, recall, F1Measure, FNList = ph.score_detector(detections, predictions)#Score the new results - Very simlar results...precision, recall, F1Measure, FNList = ph.score_detector(detections, predictionsNew)###Part 3.3 - So far we've only tested on the same image we are using for training, that is biased. It's like using the real exam for a practice test!#Let's test a novel image###Test the model on a new imageimgB, detectionsB, labelsB = ph.load_data_set(data_path,"img84")#This image looks very different. Can you spot the cells? Don't cheat :Dpyplot.imshow(imgB)if(VIEW):    pyplot.show()#Let's see the answerpyplot.imshow(imgB)pyplot.autoscale(False)pyplot.plot(detectionsB[:, 0], detectionsB[:, 1], 'kx')if(VIEW):    pyplot.show()#TODO 3.3- Get the features, run the prediction, and score the result for the new image for both#classif_oneTemplate and classif_multipleTemplates. Follow the examples above.#How do the results compare?#This image looks very different and so our templates didn't match very well ph.display_detection_results(imgB, detectionsB, predictionsB, FNList, VIEW)###Part 3.4 - We are on the right track, but we want to simulate your medical training. You don't just look at one image, you see many!###Setup the same experiment but train on the first N images with K templates per image#I suggest the following values:max_template_count = 10img_count = 23N = 3500 * img_counttemplate_img_count = 15#First get all of the templatestemplate_set = ph.get_template_set(data_path, template_img_count, max_template_count)#Get all of the features for the training imagesfeatures_train,labels_train = ph.get_features_for_image_set(data_path, img_count, N, template_set)#Visualize our new templatesph.display_templates(template_set, VIEW)#This will take a while, time for a coffeeclassif,scaler = ph.train_model(features_train, labels_train)#TODO 3.4 - Get the features, run the prediction, and score the result for the new image#Score and display the resultsprint("New results")precision, recall, F1Measure, FNList = ph.score_detector(detectionsB, predictionsBNew)ph.display_detection_results(imgB, detectionsB, predictionsBNew, FNList, VIEW)