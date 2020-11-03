
from django.shortcuts import render
from django.http import HttpResponse
from .models import Options, EndPoint, Ttest, ChiSQ , Anova , LogRank, Anova2
from django.template import RequestContext
from scipy.stats import nct,ncx2,chi2, ncf
from scipy.stats import t , f
from numpy import *
import numpy as np
import math

# Create your views here.

 
def home(request):
    return render(request,"home.html");


def add(request):
        val1= int(request.GET['both']);
        return render(request,"result.html",{'result':val1});
        
        
def step2(request):
       
            val2= int(request.GET['endeffect']);
            if val2==0:
               return render(request,"result.html",{'result':val2});
            elif val2==1:
               return render(request,"step9.html",{'result2':val2});
            elif val2==2:
               return render(request,"dichot.html",{'result2':val2});
            elif val2==3:
               return render(request,"survival.html",{'result2':val2});
            elif val2==4:
                 return render(request,"count.html",{'result2':val2});
            elif val2==5:
                 return render(request,"ttest.html",{'result2':val2});
            elif val2==6:
                 return render(request,"multicat.html",{'result2':val2});
            elif val2 == 7:
                 return render(request,"measure.html",{'result2':val2});
#******************************************************************************************************************            
#******************************************************************************************************************        
def EndpointsWithSex(request):
       
            val2= int(request.GET['sexeffect']);
            if val2==1:
               return render(request,"Sexdichot.html",{'result2':val2});
            elif val2==2:
               return render(request,"Sexsurvival.html",{'result2':val2});
            elif val2==3:
                 return render(request,"Sexmeasurettest.html",{'result2':val2});
            elif val2==5:
                 return render(request,"Sexmulticat.html",{'result2':val2});
            elif val2 ==4:
                 return render(request,"Sexmeasure2Anova.html",{'result2':val2});    
#******************************************************************************************************************            
#---------------------------------------------------------------------------------------------
#   This function calculates the power of a two-sample, two-tailed
#   t-test given:
#
#   Significance level alpha
#   Assumed population means mu_1 and mu_2
#   Assumed population standard deviations sigma_1 and sigma_2
#   Sample sizes n_1 and n_2
#---------------------------------------------------------------------------------------------

def t_test_power(alpha,mu_1,mu_2,n_1,n_2,sigma_1,sigma_2):
    nu = ((sigma_1**2/n_1+sigma_2**2/n_2)**2)/(((sigma_1**2/n_1)**2)/(n_1-1)+((sigma_2**2/n_2)**2)/(n_2-1))

    theta = (mu_1-mu_2)/((sigma_1**2/n_1+sigma_2**2/n_2)**0.5)

    t_cut = t.ppf(1-alpha/2,nu)
    
    power = 1-(nct.cdf(t_cut,nu,theta)-nct.cdf(-t_cut,nu,theta))
    
    return power
#---------------------------------------------------------------------------------------------
#---This function uses t_test_power to determine the required sample size
#   by simply starting at a low n and incrementing.  Right now, it assumes
#   n_1=n_2
#
#---If alpha, power, sigma_1 or sigma_2 are out of range, the function
#   returns 0, which is an error code
#---------------------------------------------------------------------------------------------

def sample_sizet(alpha,mu_1,mu_2,sigma_1,sigma_2,power):

    n = 0
    the_power = 0
    if alpha>0 and alpha<1 and power>0 and power<1 and sigma_1>0 and sigma_2>0:
        n = 3
        while the_power<power and n<10000:
            n = n+1
            n_1 = n
            n_2 = n
            the_power = t_test_power(alpha,mu_1,mu_2,n_1,n_2,sigma_1,sigma_2)
    return n
#---------------------------------------------------------------------------------------------
# This functions calculates the inputs for the t-test.
#---------------------------------------------------------------------------------------------    
def ttest(request):

        ttest1= Ttest()
        
        ttest1.alpha= float(request.GET['alpha']);
        ttest1.arm1 = request.GET['arm1'];
        ttest1.arm2 = request.GET['arm2'];
        ttest1.mu_1= float(request.GET['mu_1']);
        ttest1.mu_2= float(request.GET['mu_2']);
        ttest1.sigma1= float(request.GET['sigma1']);
        ttest1.sigma2= float(request.GET['sigma2']);
        ttest1.power= float(request.GET['power']);
        
        ttest1.n_per_arm = sample_sizet(ttest1.alpha,ttest1.mu_1,ttest1.mu_2,ttest1.sigma1,ttest1.sigma2,ttest1.power);
           
        return render(request,"ttest_res.html",{'ttest1':ttest1});
#******************************************************************************************************************                        
#******************************************************************************************************************            
#******************************************************************************************************************            
#   This function calculates the power of a chi-squared test for the
#	independence of the rows and columns of an RxC contingency table
#	given
#
#   Significance level alpha
#   Alternative hypothesis relative frequencies pi_1
#   Total Sample size N
#   power is returned
#---------------------------------------------------------------------------------------------
def chisq_test_power(alpha,pi_1,N):

	R = pi_1.shape[0]
	C = pi_1.shape[1]
	pi_0 = np.zeros( (R, C) )
	w = 0
	for r in range(0,R):
		for c in range(0,C):
			pi_0[r][c] = np.sum(pi_1[r])*np.sum(pi_1[:,c])
			w = w + (pi_1[r][c]-pi_0[r][c])*(pi_1[r][c]-pi_0[r][c])/pi_0[r][c]

	w = w*N

	chisq_cut = chi2.ppf(1-alpha,(R-1)*(C-1))
	
	power = 1-ncx2.cdf(chisq_cut,(R-1)*(C-1),w)
	
	return power
#---------------------------------------------------------------------------------------------
#---This function uses chisq_test_power to determine the required sample size
#   by simply starting at a low N and incrementing. 
#
#	alpha significance level of test
#	pi_1 array of probabilities undel alternative hypothesis
#	power desired power
#
#---If alpha or power are out of range, or pi_1 does not sum to 1, the function
#   returns 0, which is an error code.  Otherwise, it returns the sample size
#---------------------------------------------------------------------------------------------
def sample_size_chisq(alpha,pi_1,power):

	N = 0
	if alpha>0 and alpha<1 and power>0 and power<1 and np.sum(pi_1)==1:
		the_power = 0
		R = pi_1.shape[0]
		C = pi_1.shape[1]
		N = R*C*2
		while the_power<power and N<10000:
			N = N+1
			the_power = chisq_test_power(alpha,pi_1,N)
	return N
#---------------------------------------------------------------------------------------------
# This functions calculates the inputs for the chisq-test.
#---------------------------------------------------------------------------------------------    
def chisq(request):

        chisq1= ChiSQ();
        
        flag = request.GET['flag'];
        
        if flag == '1':
            chisq1.matrix1 = request.GET['valuelist1'];
            chisq1.matrix2 = request.GET['valuelist2'];
            entries1 = list(map(float,chisq1.matrix1.split()));
            entries2 = list(map(float,chisq1.matrix2.split()));
            entries = entries1 + entries2;
        else:    
            chisq1.matrix = request.GET['valuelist'];
            entries = list(map(float,chisq1.matrix.split()));
        
        chisq1.alpha= float(request.GET['alpha']);
        chisq1.power= float(request.GET['power']);
        
        chisq1.rows=int(request.GET['rows']);
        chisq1.cols=int(request.GET['cols']);
        
        pi_1= np.array(entries).reshape(chisq1.rows,chisq1.cols);
#        pi_1= np.array(entries).reshape(4,4);
        
        chisq1.pi_1 = pi_1/np.sum(pi_1);
        chisq1.n_per_arm = sample_size_chisq(chisq1.alpha,chisq1.pi_1,chisq1.power);
          
        return render(request,"multicat_res.html",{'chisq1':chisq1});
#---------------------------------------------------------------------------------------------
#******************************************************************************************************************            
#******************************************************************************************************************            
#******************************************************************************************************************            
#---------------------------------------------------------------------------------------------
#   This function calculates the power of an f-test given:
#
#   Significance level alpha
#   Assumed vector of population means mus
#   Assumed population standard deviation sigma
#   Common sample size n/arm n

#---------------------------------------------------------------------------------------------
def anovaf_test_power(arms,alpha,mus,n,sigma):
    ndf = arms-1
    
    ddf = (arms-1)*(n)

    nu = n*sum((mus-np.mean(mus))**2)/sigma/sigma
    
    f_cut = f.ppf(1-alpha,ndf,ddf)
    
    power = 1-ncf.cdf(f_cut,ndf,ddf,nu)
    
    return power
#---------------------------------------------------------------------------------------------
#   This function uses f_test_power to determine the required sample size
#   by starting at a low n (3) and incrementing.  
#
#   If alpha, power or sigma are out of range, the function
#   returns 0, which is an error code
#---------------------------------------------------------------------------------------------
def sample_size_anova(arms,alpha,mus,sigma,power):
    n = 0
    if alpha>0 and alpha<1 and power>0 and power<1 and sigma>0:
        the_power = 0
        n = 2
        while the_power < power and n<10000:
            n = n+1
            the_power = anovaf_test_power(arms,alpha,mus,n,sigma)
            
    return n
#---------------------------------------------------------------------------------------------
# This functions calculates the inputs for the anova test
#---------------------------------------------------------------------------------------------  
def anova(request):

        anova1= Anova();
        
        entries=[];
        flag = request.GET['flag'];
        
        if flag == '1':
            anova1.matrixf = request.GET['valuelistf']
            entries = list(map(float,anova1.matrixf.split()))
            
        anova1.alpha= float(request.GET['alpha']);
        anova1.matrix = request.GET['valuelist'];
        anova1.power=float(request.GET['power']);
        anova1.sigma=float(request.GET['sigma']);
        anova1.arms=float(request.GET['arms']);
		
        entries= entries + (list(map(float,anova1.matrix.split())));
        anova1.mus = np.array([entries],dtype=float);
        anova1.len=len(anova1.mus);
        anova1.n_per_arm = sample_size_anova( anova1.arms,anova1.alpha,anova1.mus,anova1.sigma,anova1.power);
        
        return render(request,"measure_res.html",{'anova1':anova1,'flag':flag});
        
#---------------------------------------------------------------------------------------------
#******************************************************************************************************************            
#******************************************************************************************************************            
#******************************************************************************************************************            
#---------------------------------------------------------------------------------------------
#   This function uses logrank_test_power to determine the required sample size
#   by simply starting at a low N and incrementing. 
#
#       alpha significance level of test
#       pi_1 array of probabilities undel alternative hypothesis
#       power desired power
#       T end of observation period
#
#   If alpha or power are out of range, the function returns 0, which is an 
#       error code.  Otherwise, it returns the sample size
#---------------------------------------------------------------------------------------------
def sample_size_log(arms,alpha,meds,power,T):

        N = 0
        if alpha>0 and alpha<1 and power>0 and power<1 and T>0.0:
                the_power = 0
                N = 3
                while the_power<power and N<10000:
                        N = N+1
                        the_power = logrank_test_power(arms,alpha,meds,N,T)
        return N

#---------------------------------------------------------------------------------------------
#---This function calculates the power of the log-rank test
#       for the equality of k survival functions 
#
#   Significance level alpha
#   Alternative median survivals meds
#   Total Sample size N
#   T end of observation period
#   power is returned
#
#---------------------------------------------------------------------------------------------
def logrank_test_power(arms,alpha,meds,N,T):

        K = arms
       
        hs = math.log(2)/meds
        lhs = np.log(hs)
        psi_2 = 0
        for k in range(0,K):
                psi_2 = psi_2+(N*(1-math.exp(-hs[0,k]*T))*(lhs[0,k]-np.mean(lhs))**2)
        chisq_cut = chi2.ppf(1-alpha,K-1)
        power = 1-ncx2.cdf(chisq_cut,K-1,psi_2)
        return power
#---------------------------------------------------------------------------------------------
def logrank(request):

        logrank= LogRank();
        flag = request.GET['flag'];
        entries=[];
        if flag == '1':
            logrank.matrixf = request.GET['valuelistf']
            entries = list(map(float,logrank.matrixf.split()))
            logrank.arms=int(request.GET['arms']);
        
        logrank.matrix = request.GET['valuelist'];
        logrank.alpha= float(request.GET['alpha']);
        logrank.power= float(request.GET['power']);
        
        logrank.observation= float(request.GET['observation']);
        
        entries = list(map(float,logrank.matrix.split())) + entries;
        logrank.meds = np.array([entries],dtype=float);
        
        if flag == '1':
            logrank.meds=reshape(logrank.meds,(2,logrank.arms));
           
        
        if flag == '0':
           logrank.arms= logrank.meds.size
        
        logrank.n_per_arm = sample_size_log(logrank.arms,logrank.alpha,logrank.meds,logrank.power,logrank.observation);
           
        return render(request,"survival_res.html",{'logrank':logrank, 'flag':flag});
#---------------------------------------------------------------------------------------------
#******************************************************************************************************************            
#******************************************************************************************************************            
#******************************************************************************************************************            
#---------------------------------------------------------------------------------------------   
#   This function calculates the power of an f-test given:
#
#   Significance level alpha
#   Matrix of assumed population means mus
#   Assumed population standard deviation sigma
#   Common sample size n/cell, n
#---------------------------------------------------------------------------------------------   
def anova2_test_power(arms,alpha,mus,n,sigma):
    n_cell=0
    ndf=0
    ddf=0
    n_cell = arms*2
    ndf = n_cell-1
    ddf = n*n_cell-ndf-1
    nu = n*np.sum((mus-np.mean(mus))**2)/sigma/sigma
    f_cut = f.ppf(1-alpha,ndf,ddf)
    power = 1-ncf.cdf(f_cut,ndf,ddf,nu)
    return power

#---------------------------------------------------------------------------------------------   
#   This function uses f_test_power to determine the required sample size
#   by starting at a low n (3) and incrementing.  
#
#   If alpha, power or sigma are out of range, the function
#   returns 0, which is an error code
#---------------------------------------------------------------------------------------------   
def sample_size_anova2(arms,alpha,mus,sigma,power):

    n = 0
    if (alpha>0.0) and (alpha<1.0) and (power>0.0) and (power<1.0) and (sigma>0.0):
        the_power = 0.0
        n = 2
        while (np.all(the_power<power) and np.all(n<10000)):
            n = n+1
            the_power = anova2_test_power(arms,alpha,mus,n,sigma)
    return n

#---------------------------------------------------------------------------------------------
# This functions calculates the inputs for the 2 WAY anova test
#---------------------------------------------------------------------------------------------  
def anova2(request):

        anova2way= Anova2();
        
        flag = request.GET['flag'];
        entries=[];
                
       
        anova2way.matrix = request.GET['valuelist'];
        entries = list(map(float,anova2way.matrix.split()))
        anova2way.matrixf = request.GET['valuelistf']   
        anova2way.alpha= float(request.GET['alpha']);
        
        anova2way.power=float(request.GET['power']);
        anova2way.sigma=float(request.GET['sigma']);
        anova2way.arms=float(request.GET['arms']);
		
        entries= entries + (list(map(float,anova2way.matrixf.split())));
        anova2way.mus = np.array([entries],dtype=float);
        
        anova2way.n_per_arm = sample_size_anova2( anova2way.arms,anova2way.alpha,anova2way.mus,anova2way.sigma,anova2way.power);
        
#        n_cell=0
#        ndf=0
#        ddf=0
#        n_cell = 6
#        ndf = n_cell-1
#        ddf = 3*(n_cell-ndf-1)
#         nu = 3*np.sum((anova2way.mus-np.mean(anova2way.mus))**2)/anova2way.sigma/anova2way.sigma
#        f_cut = f.ppf(1-anova2way.alpha,ndf,ddf)
#        anova2way.test=f_cut
#        power = 1-ncf.cdf(f_cut,ndf,ddf,nu)
       
        
        return render(request,"measure_res.html",{'anova2way':anova2way,'flag':flag});