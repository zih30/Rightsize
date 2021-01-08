from django.db import models


# Create your models here.
class Ttest:
    alpha:float
    arm1:str
    arm2:str
    mu_1:  float
    mu_2:  float
    sigma1: float
    sigma2: float
    power:float
    n_per_arm:int
    
class ChiSQ:
    alpha:float
    power:float
    n_per_arm:int
    matrix:float
    matrix1:float
    matrix2:float
    rows:int
    cols:int
    
class LogRank:
    alpha:float
    power:float
    n_per_arm:int
    matrix:str
    arms:int

    
class Anova:
    alpha:float
    power:float
    n_per_arm:float
    sigma:float
    matrix:str
    arms:int
    
class Anova2:
    alpha:float
    power:float
    n_per_arm:float
    sigma:float
    matrix:str
    matrixf:str
    arms:int
    test:float
    
    
class Options:
     yes:bool
     no:bool

class EndPoint:
    Dichot:bool
    Survival:bool
    Count:bool
    Measure:bool
    MultiCat:bool
    
