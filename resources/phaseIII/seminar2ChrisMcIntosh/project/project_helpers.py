import numpy as npfrom scipy import ndimagefrom io import StringIOfrom skimage import iofrom skimage.feature import match_templatefrom skimage import colorfrom scipy.spatial import distance as dtfrom matplotlib import pyplotfrom skimage import filtersfrom skimage.feature import peak_local_maxfrom skimage.color import rgb2hedfrom sklearn.ensemble import RandomForestClassifierfrom joblib import Parallel, delayedimport multiprocessingfrom multiprocessing import Pooldef load_data_set(data_path, image_name):    """Loads a data set."""    img = io.imread(data_path + "Detection/" + image_name + "/" + image_name + ".bmp")        tmp = np.loadtxt(data_path + "Detection/" + image_name + "/" + image_name + "_detection.txt",ndmin=2)    detections = tmp[:, 0:2] - 1    class_labels = tmp[:, 2]    return [img, detections, class_labels]def get_training_samples(featuresND, detections, img_gray, N):    """Get training samples from a set of features. One sample per detection,    and a bunch of random samples where there is no nuclei centre."""    N = int(N)    #Positive examples are detections, negative examples are other pixels    samples = detections.shape[0]    num_features = featuresND.shape[2]    features_train = np.zeros([N, num_features])    labels_train = np.zeros([N, 1])    #Negative examples are where there is not a detection    negative_sample = np.ones(img_gray.shape)    k = 3    for idx in range(0, samples):        features_train[idx,:] = featuresND[int(detections[idx, 1]), int(detections[idx, 0]),:]        labels_train[idx] = 1        #Mark a 7 pixel by 7 pixel mask in cell around centroid        yC = int(detections[idx, 1])        xC = int(detections[idx, 0])        #Account for boundary conditions where the mask extends beyond the image        miny = max((0, yC-k))        maxy = min((img_gray.shape[1]-1, yC+k+1))        minx = max((0, xC-k))        maxx = min((img_gray.shape[0]-1, xC+k+1))        negative_sample[miny:maxy, minx:maxx] = 0     #Get all of the negative sample locations    NS = np.where(negative_sample == 1)    count = 0    for fx in np.random.choice(NS[0].shape[0], N-samples,replace=False):        idy = NS[0][fx]        idx = NS[1][fx]        if ((negative_sample[idy, idx] == 1)):            features_train[samples + count] = featuresND[idy, idx, :]            labels_train[samples + count] = -1            count = count + 1    return features_train, labels_traindef process_input(idx, img_gray, templates):        # TODO 3.1.5 Fill in the return value of the process_input function to match        # the idx'th template in templates using the template matching from skikit-image.        # store the result in template_Response. Set optional input pad_input=True            #template_Response = match_template(img_gray, templates[:,:,idx], pad_input=True)    return template_Responsedef calculate_features(img, templates):    "This function calculates the image features, 1 per voxel"    # Note, this code will run as-is. However, as you add more    # and more features it will increase the accuracy.        # For the rest of the features first convert the image to a    # floating point array (it is currently integer valued).    img = img.astype('float')    #When finished implementing 3.1, remove this line    featuresND = img        #TODO PART 3.1.x fill in the missing code using the imageProcessing_Tutorial as a guide        #TODO 3.1.1 - Replace rgb2gray below with a conversion to HED stain using skikit-image. Extract the Haematoxylin channel, and    #call the result img_gray, it should be MxN    img_gray = my_rgb2gray(img)    #TODO 3.1.1 thru 3.1.4 - You are given a simple averaging filter, F, fill in the features below    #Use the convolution function from scipy to perform convolution and calcluate the features    # Set optional parameters: mode='constant', cval=0.0 to handle the boundary conditions    F = (1 / 9) * np.array([[1, 1, 1],[1, 1, 1],[1, 1, 1]])        #Feature A - Compute the average intensity of the resulting single-channel Haematoxylin image    featureA = ndimage.convolve(img_gray, F, mode='constant', cval=0.0)    #TODO 3.1.2    #Feature B - The average Red value    featureB = ndimage.convolve(img[:,:,0].astype('float'), F, mode='constant', cval=0.0)    #TODO 3.1.3    #Feature C - The average Green value    featureC = ndimage.convolve(img[:,:,1].astype('float'), F, mode='constant', cval=0.0)    #TODO 3.1.4    #Feature D - The average Blue value    featureD = ndimage.convolve(img[:,:,2].astype('float'), F, mode='constant', cval=0.0)    #TODO 3.1 - Stack all of the features into a MxNxF array where F is the number of features for an MxN image    featuresND = np.dstack((featureA, featureB, featureC, featureD))    #Ensure that a single template is still a 3-dimensional array    if templates.ndim <= 2:        templates = np.reshape(templates, (templates.shape[0], templates.shape[1], 1))    template_features = np.zeros((featuresND.shape[0], featuresND.shape[1], templates.shape[2]))    num_cores = multiprocessing.cpu_count()         #Parallel(n_jobs=2,backend="threading")(delayed(test)(i ** 2) for i in range(100))        # This code calls process_input to compute the template response across all templates, you do not need to modify this    template_features = Parallel(n_jobs=num_cores,backend="threading")(delayed(process_input)(idx, img_gray, templates) for idx in range(0, templates.shape[2]))    template_features = np.asarray(template_features)    template_features = np.rollaxis(template_features,0,3)    # Stack all of the features for processing        featuresND = np.dstack((featuresND, template_features))    return featuresNDdef get_templates(img_gray, detections, max_template_count):    k = 10    templates = np.zeros([2*k+1, 2*k+1, max_template_count])    count = 0    #Get templates for a random set of detections    for idx in np.random.choice(detections.shape[0], min(detections.shape[0], max_template_count),replace=False):        ###################        #TODO 3.2 (BELOW) - Get each template, centred at the detection pixel        # You will get a box (array) of dimensions 2*k+1, 2*k+1 centred at the pixel        #Begin by drawing by hand an example pixel at [20,30] and a box of width 21, with 10 pixels to the        #left and 10 to the right        #Psudecode:        #1) Get the pixel location from detections in x,y co-ordinate and convert it to integer.        #1b) Remember, the image is accessed as img_gray[yValue,xValue]        #2) Ensure that the box will fit in the image, for example, that the right edge of the box is not > img_gray.shape[0]        #2b) There are a total of 4 such boundary conditions to ensure the box is within the left, right, bottom, and top boundaries of the image.        #3) Only if the box is within the image, access the img_gray array and retreive the 2*k+1, 2*k+1 elements        #4) Store them in templates along the 3rd dimension        #5) Increment a counter so that you're filling up the templates array from 0 to end, along its 3rd dimension        #If you get stuck, draw the example using the convolution example from the lecture slides        ###################          xC = int(detections[idx,0])        yC = int(detections[idx,1])        if( (xC > k) and ((xC+k+1) < img_gray.shape[1]) and (yC > k) and ((yC+k+1) < img_gray.shape[0])):            T = img_gray[yC-k : yC+k+1, xC-k : xC+k+1]            templates[:,:,count] = T            count = count +1    #TODO 3.2: remove any un-used template slots from the end of your array down to your counter    templates = templates[:,:,0:count]    return templatesdef predict_centres(img, features_img, scaler,classif):    t = features_img.reshape([img.shape[0] * img.shape[1], features_img.shape[2]])    t = scaler.transform(t)    probability_prediction = classif.predict_proba(t)[:,1]        probability_prediction = probability_prediction.reshape([img.shape[0], img.shape[1]])    probability_prediction = filters.gaussian(probability_prediction, sigma=0.6)    probability_prediction = probability_prediction / np.max(probability_prediction)    centres = peak_local_max(probability_prediction, min_distance=6, threshold_abs=0.3)    centres = np.fliplr(centres)    return probability_prediction, centresdef score_detector(detections_ground_truth, detections_proposed):    """Find the number of true and false positive/negatives, and score them."""    #Set the threshold to within 12 pixels    threshold = 12    P = detections_ground_truth.shape[0]    TP = 0    FN = 0    FNList = np.zeros((P,1))     #TODO Part 2.1    #Fill in the algorithm to count true positives and false negatives    #For each ground truth detection, see if there is a proposed detection within a threshold distance    #Keep track of the indexes of detectionGroundTruth that are false negative (i.e. missed) in FNList    #Psuedocode    #For each ground truth point    #   Get the Euclidean distance between it and all detections_proposed points    #   Get the minimum distance    #   If the distance is less than threshold, it is a true positve; Remove that point from detections_proposed since it can only be true once    #   Else, it is a false negative        #Get the distance between all pairs of points    dst = dt.cdist(detections_ground_truth, detections_proposed, 'euclidean')    for idx in range(0,P):        min_loc = np.argmin(dst[idx,:])        if ((dst[idx,min_loc]) <= threshold):            TP = TP + 1            dst[:,min_loc] = 2*threshold # Ensure it cannot be chosen again        else:            FN = FN + 1            FNList[idx] = 1    #Every proposed detection that isn't a TP is a FP    FP = detections_proposed.shape[0] - TP    #TODO Part 2.1 - Calculate the metrics    precision = TP / (TP + FP)    recall = TP / P    F1Measure = 2 * TP / (2 * TP + FP + FN)    print('Precision: %1.3f; Recall %1.3f; F1Measure %1.3f'%(precision,recall,F1Measure))      return precision, recall, F1Measure, FNListdef display_detection_results(img, detections, predictions, FNList, VIEW):    pyplot.imshow(img)    pyplot.autoscale(False)    pyplot.plot(predictions[:, 0], predictions[:, 1], 'go')    NS = np.where(FNList == 1)    pyplot.plot(detections[NS[0], 0], detections[NS[0], 1], 'rx')    NS = np.where(FNList == 0)    pyplot.plot(detections[NS[0], 0], detections[NS[0], 1], 'kx')    if(VIEW):        pyplot.show()def display_templates(templates, VIEW):    """Visualize all of the templates."""    N = 10    count = 0    for idx in range(0,(templates.shape[2]//(N*N))):        f, axarr = pyplot.subplots(N, N)        for idxx in range(0, N):            for idxy in range(0, N):                if(count < templates.shape[2]):                    axarr[idxx,idxy].imshow(templates[:,:,count], cmap=pyplot.cm.gray)                    axarr[idxx,idxy].axis('off')                    count = count + 1        if(VIEW):            pyplot.show(block=False)            input("Press Enter to continue...")            pyplot.close(f)        def train_model(features_train,labels_train):    #Rescale the data to zero mean and unit variance       from sklearn import preprocessing    scaler = preprocessing.StandardScaler().fit(features_train)    features_train_scale = scaler.transform(features_train)    classif = RandomForestClassifier(random_state=0, min_samples_split=1, class_weight='balanced', n_jobs=-1)    classif.fit(features_train_scale, np.ravel(labels_train))     return classif, scaler   def get_template_set(data_path, img_count, max_template_count):    for idx in range(1, img_count + 1):        #Load an image        img_name = "img%.0f"%(idx)        img_local, detections_local, labels_local = load_data_set(data_path, img_name)        img_gray_local = my_rgb2gray(img_local)        templates_local = get_templates(img_gray_local, detections_local, max_template_count)        if(idx == 1):            template_set = templates_local        else:            template_set = np.dstack((template_set, templates_local))    return template_setdef get_features_for_image_set(data_path, img_count, N, template_set):    """Get a total of N features for images 1:img_count."""    #TODO 3.4 - Write a loop to extract the features from the images,     # get the training samples, and stack the results into features_train    # and labels_train.    for idx in range(1, img_count + 1):        #Load an image        img_name = "img%.0f"%(idx)        img_local, detections_local, labels_local = load_data_set(data_path, img_name)        img_gray_local = my_rgb2gray(img_local)        featuresLocal= calculate_features(img_local,template_set)        features_train_local, labels_train_local = get_training_samples(featuresLocal, detections_local, img_gray_local, N / img_count)        print('Finished extracted features from image', idx)        if(idx == 1):            features_train = features_train_local            labels_train = labels_train_local        else:            features_train = np.vstack((features_train, features_train_local))            labels_train = np.vstack((labels_train, labels_train_local))    return features_train, labels_traindef my_rgb2gray(img):    return color.rgb2hed(img)[:,:,0]def validate_image(imgB,detectionsB,template_set, scaler, classif):    #imgB - the image to extract features and validate    #detectionsB - the detections for the image    #template_set - the set of templates used by the model to build features    #scaler,classif - the learned classifier variables    featuresImgB = ph.calculateFeatures(imgB,template_set)    probabilityPredictionNew,predictionsBNew = ph.predict_centres(imgB,featuresImgB,scaler,classif)    precision_local,recall_local,F1Measure_local,tmp = ph.score_detector(detectionsB,predictionsBNew)    precision = precision_local    recall = recall_local    F1Measure = F1Measure_local        return precision, recall, F1Measuredef validate_model(test_set, data_path, template_set, scaler, classif):    #test_set - numpy array of image numbers to test as an Nx1 array    #data_path - the file path to the data set    #template_set - the set of templates used by the model to build features    #scaler,classif - the learned classifier variables    num_test_images = test_set.shape[0]    precision = np.zeros((num_test_images,1))    recall = np.zeros((num_test_images,1))    F1Measure = np.zeros((num_test_images,1))    count = 0    for idx in range(0,num_test_images):        img_name = "img%.0f"%(test_set[idx])        imgB, detectionsB, labelsB = ph.load_data_set(data_path,img_name)        featuresImgB = ph.calculateFeatures(imgB,template_set)        probabilityPredictionNew,predictionsBNew = ph.predict_centres(imgB,featuresImgB,scaler,classif)        precision_local,recall_local,F1Measure_local,tmp = ph.score_detector(detectionsB,predictionsBNew)        precision[count] = precision_local        recall[count] = recall_local        F1Measure[count] = F1Measure_local        count = count +1    print("Precision %f +/- %f; Recall %f +/- %f; F1Measure %f +/- %f"%(np.mean(precision),np.std(precision),np.mean(recall),np.std(recall),np.mean(F1Measure),np.std(F1Measure)))    return precision, recall, F1Measure