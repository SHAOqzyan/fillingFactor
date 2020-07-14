
import numpy as np
from scipy.optimize import leastsq




########
def ffFunction(x,a, b, c):
    #return a*(x-x0)**2+c
    #print a, b,c
    return a*np.exp(-b*x)  +c

def dya(x,a,b,c):
    return np.exp(-b*x)

def dyb(x,a,b,c):
    return a*(-b)*np.exp(-b*x)

def dyc(x,a,b,c):
    return 10






class myBffLsq(object):

    def __init__(self):
        pass


    def residual(self,variables,x,data, dataError):
        """

        :param variables:
        :param x:
        :param data:
        :param dataError:
        :return:
        """
        #a,b,c=variables
        model=ffFunction(x,*variables)
        #return   (model- data )/ dataError
        return   ( data -model)/ dataError

    def Jacobi(self,variables,x,data,dataError):

        """

        :param variables:
        :param x:
        :param data:
        :param dataError:
        :return:
        """

        #the jacobi matrix is N*3, because we have three parameters

        partialA=  dya(x,*variables)/dataError
        partialB=  dyb(x,*variables)/dataError
        partialC=  dyc(x,*variables)/dataError



        Jr=np.asarray( [ partialA,  partialB  , partialC    ] )


        return Jr



    def doLSQ(self,x,y,yError):

        """

        :param x:
        :param y:
        :param yError:
        :return:
        """

        variables=[np.mean(y), 0.5,np.mean(y)]

        out=leastsq( self.residual,variables,args=(x,y,yError) )
        print out ,"????????????????"

        #can you calculate the