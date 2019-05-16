import numpy as np
import json
from sklearn.externals import joblib




def predict_los(x1,x2,x3,x4,x5):
    X=[[x1,x2,x3,x4,x5]]
    try:
        model=joblib.load('mymodel.model')
    except FileNotFoundError:
        print("模型不存在")
    Y=model.predict(X)#验证
    return Y[0]
