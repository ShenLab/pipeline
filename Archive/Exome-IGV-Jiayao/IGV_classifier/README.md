Scripts for convolutional nerural network IGV image classifier

the classifier takes IGV snapshot image as input and predict the image into 5 classes: snp, insertion, deletion, complex, failed IGV and uncertain.
The accuracy is higher among snps and there still have many false positive among indels(~15%)


### Requirement
[PIL] (https://github.com/python-pillow/Pillow)

[keras] (https://github.com/fchollet/keras)


### Tutorials

- [preprocess one image] (process_images.ipynb)

- [prepare for training data] (prepare_training_data.ipynb)

- [train a convolution neural network model] (IGV_cnv_train.ipynb)

- [IGV prediction](IGV_prediction.ipynb)
