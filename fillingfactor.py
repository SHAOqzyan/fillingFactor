
import numpy as np
import matplotlib as mpl
mpl.use('agg')
import scipy.odr.odrpack as odrpack
import radio_beam
from spectral_cube import SpectralCube
from astropy import units as u

from mwispDBSCAN import MWISPDBSCAN

import matplotlib.pyplot as plt
from myPYTHON import *
from astropy.io import fits
import glob
import os
import sys
import seaborn as sns
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from scipy import optimize
from progressbar import *
from astropy.table import Table,vstack
import gc

#sys.path.insert(1, '/home/qzyan/WORK/myDownloads/MWISPcloud')
#sys.path.insert(1, '/home/yqzpmo/github/MWISPdbscan')

#from onlyDBSCAN import allDBSCAN

#doAllDBSCAN = allDBSCAN()

doFITS=myFITS() #used to deal fits with myPYTHON


########
def ffFunction(x,a, b, c):
    #return a*(x-x0)**2+c
    #print a, b,c
    return a*np.exp(-b*x)  +c


def ffFunctionOdr(B,x ):
    a,b,c=B

    #return a*(x-x0)**2+c
    return a*np.exp(-b*x)  +c




##########################################

def dya(x,a,b,c):
    return np.exp(-b*x)

def dyb(x,a,b,c):
    return a*(-b)*np.exp(-b*x)

def dyc(x,a,b,c):
    return 1
####################################

def ffAndSizeFunc(B, x ):
    #return  b- b*np.exp(-a*x )
    #return  1- b*np.exp(-a*x)

    #return  B[1]- B[1]*np.exp(-B[0]*x)
    #return 1- 1*np.exp(-B[0]*x)
    a,b=B
    return  a-b/x

    #return   np.exp(-a/np.sqrt(x) )


def testffAndSizeFunc(x,  a  ):
    #return  b- b*np.exp(-a*x )
    #return  1- b*np.exp(-a*x)

    #return  b- b*np.exp(-a*x )
    #return  a-b/(x-0.564)
    #return  a-b/(x+c)
    zeroPoint= 0 # 0.564  #0.977

    #c=b/a  -zeroPoint


    #c=b
    return x/(x +a)


def testffAndSizeFunc1(x, a):
    return x / (x + a)

def testffAndSizeFunc2(x, a,b):
    return a- a*np.exp(-b*x )

def testffAndSizeFunc3(x, a,b,c):
    return a*np.exp(-b*x)  +c


    #return   np.exp(-a/np.sqrt(x) )
def testffAndSizeFuncODR(B,x  ):
    #return  b- b*np.exp(-a*x )
    #return  1- b*np.exp(-a*x)

    #return  b- b*np.exp(-a*x )
    #return  a-b/(x-0.564)
    #return  a-b/(x+c)

    return x/(x +B[0])

    #return   np.exp(-a/np.sqrt(x) )



class checkFillingFactor(object):

    rootPath="./"

    #get the curr


    saveFITSPath="/home/qzyan/WORK/diskMWISP/fillingFactorData/"
    #saveFITSPath="./"

    dataPath= saveFITSPath+"data/"
    rawFITS=  dataPath + "G2650Local30.fits"
    rmsFITSPath= rootPath+"rmsFITSpath/"
    #tmpPath= rootPath +"tmpFiles/"
    tmpPath= saveFITSPath +  "tmpFiles/"

    figurePath=  rootPath +"figurePath/"
    statPath=rootPath+"statPath/"
    intFigurePath = rootPath+"cloudIntMap/"

    paperFigurePath=rootPath+"paperFigurePath/"

    ffFigurePath=rootPath+"ffFigures/"


    cloudCubePath= saveFITSPath+ "cloudCubes/"

    #badChannel spectral path

    badChannelSpecPath=rootPath+"badChannelSpec/"
    #maskTmpPath=saveFITSPath+"maskedFITS/"

    ######## Out
    codeOutCO12="OutCO12"
    outCO12FITS =dataPath+ "OutCO12.fits"

    codeOutCO13="OutCO13"
    outCO13FITS =dataPath+ "OutCO13.fits"

    codeOutCO18="OutCO18"
    outCO18FITS =dataPath+ "OutCO18.fits"

    #raw
    codeRawOutCO12="rawOutCO12"
    rawOutCO12FITS =dataPath+ "rawOutCO12.fits"

    codeRawOutCO13="rawOutCO13"
    rawOutCO13FITS =dataPath+ "rawOutCO13.fits"

    codeRawOutCO18="rawOutCO18"
    rawOutCO18FITS =dataPath+ "rawOutCO18.fits"




    ############## Local
    codeLocalCO12="LocalCO12"
    localCO12FITS =dataPath+ "LocalCO12.fits"


    codeLocalCO13="LocalCO13"
    localCO13FITS =dataPath+ "LocalCO13.fits"


    codeLocalCO18="LocalCO18"
    localCO18FITS =dataPath+ "LocalCO18.fits"

    ############## rawLocal
    codeRawLocalCO12="rawLocalCO12"
    rawLocalCO12FITS =dataPath+ "rawLocalCO12.fits"


    codeRawLocalCO13="rawLocalCO13"
    rawLocalCO13FITS =dataPath+ "rawLocalCO13.fits"


    codeRawLocalCO18="rawLocalCO18"
    rawLocalCO18FITS =dataPath+ "rawLocalCO18.fits"



    ########### Sgr

    codeSgrCO12= "sgrCO12"
    sgrCO12FITS =dataPath+ "sgrCO12.fits"

    codeSgrCO13 = "sgrCO13"
    sgrCO13FITS = dataPath+ "sgrCO13.fits"

    codeSgrCO18 = "sgrCO18"
    sgrCO18FITS = dataPath+ "sgrCO18.fits"

    #raw
    codeRawSgrCO12= "rawSgrCO12"
    rawSgrCO12FITS =dataPath+ "rawSgrCO12.fits"

    codeRawSgrCO13 = "rawSgrCO13"
    rawSgrCO13FITS = dataPath+ "rawSgrCO13.fits"

    codeRawSgrCO18 = "rawSgrCO18"
    rawSgrCO18FITS = dataPath+ "rawSgrCO18.fits"





    ###### Scu

    codeScuCO12= "scuCO12"
    scuCO12FITS = dataPath+ "scuCO12.fits"

    codeScuCO13 = "scuCO13"
    scuCO13FITS = dataPath+ "scuCO13.fits"

    codeScuCO18 = "scuCO18"
    scuCO18FITS = dataPath+ "scuCO18.fits"

    #raw
    codeRawScuCO12= "rawScuCO12"
    rawScuCO12FITS = dataPath+ "rawScuCO12.fits"

    codeRawScuCO13 = "rawScuCO13"
    rawScuCO13FITS = dataPath+ "rawScuCO13.fits"

    codeRawScuCO18 = "rawScuCO18"
    rawScuCO18FITS = dataPath+ "rawScuCO18.fits"



    #raw

    #not exactly
    outVrange= [-79,-6] #km/s
    localVrange= [-6,30] #km/s
    sgrVrange= [30,70] #km/s
    scuVrange= [70,120] #km/s



    #tmpPath=  "/home/qzyan/WORK/dataDisk/fillingFactorData/"

    #rawRMS = 0.5
    MWISPrmsCO12 = 0.423  ##0.49  ##K
    MWISPrmsCO13 = 0.228  ## K
    MWISPrmsCO18 = 0.229 # #K

    MWISPrmsRawCO12 = 0.485  ##0.49  ##K
    MWISPrmsRawCO13 = 0.262  ## K
    MWISPrmsRawCO18 = 0.263 # #K



    #smoothFactors =  np.arange(1.0,15.5,0.5  ) #compared to beam size, this is also the distance factor, that push cloud away by this factor
    smoothFactors =  np.arange(1.0, 10.5, 0.5  ) #compared to beam size, this is also the distance factor, that push cloud away by this factor


    noiseFactors =  np.arange(0.0, 2.1, 0.1  ) #  in Kelvin
    cutoffFactors =  np.arange(2.0,10.5,0.5  ) #  in Kelvin


    #just for test
    #smoothFactors = np.arange(1.0, 2,   0.5)  # compared to beam size, this is also the distance factor, that push cloud away by this factor
    #noiseFactors = np.arange(1.0, 2, 0.5)  # compared to  raw data rms

    calCode=""
    noiseStr="Noise"

    rawBeamSizeCO12 =   49. / 60
    rawBeamSizeCO13 =   52. / 60
    rawBeamSizeCO18 =   52. / 60



    normNoise = mpl.colors.Normalize(vmin= min(noiseFactors) , vmax= max(noiseFactors) )
    noiseColor = plt.cm.ScalarMappable(norm=normNoise, cmap=plt.cm.jet)

    normSmooth = mpl.colors.Normalize(vmin= min(smoothFactors) *  rawBeamSizeCO12 , vmax= max(smoothFactors)*rawBeamSizeCO12 )
    smoothColor = plt.cm.ScalarMappable(norm=normNoise, cmap=plt.cm.jet)


    #####################

    noiseList= [ MWISPrmsCO12, MWISPrmsCO13,   MWISPrmsCO18  ]

    #armStrList=["local", "Sagittarius", "Scutum-Centaurus","Outer" ] #the perseus arm are incorperated into the local arm

    lineStrList=["\\cofs",  "\\coss", "\\cots" ]

    localCodeList=[   codeLocalCO12,  codeLocalCO13,   codeLocalCO18  ]
    sgrCodeList=  [   codeSgrCO12,    codeSgrCO13,     codeSgrCO18    ]
    scuCodeList=  [   codeScuCO12,    codeScuCO13,     codeScuCO18    ]
    outCodeList=  [   codeOutCO12,    codeOutCO13,     codeOutCO18    ]


    rawLocalCodeList=[   codeRawLocalCO12,  codeRawLocalCO13,   codeRawLocalCO18  ]
    rawSgrCodeList=  [   codeRawSgrCO12,    codeRawSgrCO13,     codeRawSgrCO18    ]
    rawScuCodeList=  [   codeRawScuCO12,    codeRawScuCO13,     codeRawScuCO18    ]
    rawOutCodeList=  [   codeRawOutCO12,    codeRawOutCO13,     codeRawOutCO18    ]


    allRawCodeList =  rawLocalCodeList+  rawSgrCodeList +  rawScuCodeList +  rawOutCodeList



    allCodeList =  localCodeList+  sgrCodeList +  scuCodeList +  outCodeList
    allNoiseList = noiseList + noiseList + noiseList + noiseList
    allArmList=  ["local"]*3+  ["Sagittarius"]*3+  ["Scutum-Centaurus"]*3+  ["Outer"]*3    #armStrList+armStrList+armStrList+armStrList

    allLineList= lineStrList+ lineStrList+lineStrList+lineStrList

    ##################### smooth factors

    ffMWISPCol="fillingFactorMWISP"
    ffMWISPErrorCol="fillingFactorMWISPError"

    
    ffCfACol="fillingFactorCfa"
    ffCfAErrorCol="fillingFactorErrorCfa"

    aCol="para_a"
    bCol="para_b"
    cCol="para_c"


    aErrorCol="error_a"
    bErrorCol="error_b"
    cErrorCol="error_c"

    ##################### cutoff factors

    cutFFffMWISPCol = "cutFFfillingFactorMWISP"
    cutFFffMWISPErrorCol = "cutFFfillingFactorMWISPError"

    cutFFffCfACol = "cutFFfillingFactorCfa"
    cutFFffCfAErrorCol = "cutFFfillingFactorErrorCfa"

    cutFFaCol = "cutFFpara_a"
    cutFFbCol = "cutFFpara_b"
    cutFFcCol = "cutFFpara_c"

    cutFFaErrorCol = "cutFFerror_a"
    cutFFbErrorCol = "cutFFerror_b"
    cutFFcErrorCol = "cutFFerror_c"



    #################################
    touchLedgeCol= "touchLedge"
    touchBedgeCol= "touchBedge"
    touchVedgeCol= "touchVedge"

    peakChannelCol="peakChannel"
    conseCol="conseCol"

    drawCodeArea="area"
    drawCodeFlux="flux"
    drawCodeSize="size"


    rmsCO12FITS="CO12RMS.fits"
    rmsCO13FITS="CO13RMS.fits"
    rmsCO18FITS="CO18RMS.fits"


    rmsRawCO12FITS="rawCO12RMS.fits"
    rmsRawCO13FITS="rawCO13RMS.fits"
    rmsRawCO18FITS="rawCO18RMS.fits"


    smoothTag="_SmFactor_"

    noiseTag="_NoiseAdd_"

    idCol= "_idx"


    VCol="v_cen"
    LCol="x_cen"
    BCol="y_cen"

    ##
    eqLwCol="eqLwCol" #km /s

    vrmsAvgSpectraCol="vrmsAvgSpectra" #km/s, vrms calculated with average spectrum
    vmeanAvgSpectraCol="vmeanAvgSpectra" #km/s, vrms calculated with average spectrum


    def __init__(self):
        pass



    def getNList(self,TBList):


        Nlist=[]

        for eachTB in TBList:

            Nlist.append( len(eachTB) )



        return Nlist



    def produceNoiseFITS(self, channelN, rmsFITS=None,inputRmsData=None):
        """
        based on the rmsFITS, produce a spectral with channelN,  and add this cube to smooth fits to smooth fits
        smooth fits, needto calculate rms fits, better use outArmFITS,
        :param rmsFITS:
        :param channelN:
        :return:
        """
        if inputRmsData is None:
            dataRMS,headRMS=doFITS.readFITS(rmsFITS)

        else:
            dataRMS = inputRmsData

        np.random.seed()

        Ny,Nx=dataRMS.shape
        Nz=channelN

        dataNoise=np.zeros( (Nz,Ny,Nx), dtype=np.float  )


        for j in range(Ny):
            for i in range(Nx):

                if np.isnan( dataRMS[j,i]):
                    continue

                dataNoise[:,j,i]=np.random.normal(0,dataRMS[j,i],Nz)


        #print dataNoise.shape

        #fits.writeto("noiseCube.fits", dataNoise )

        return dataNoise

    def getRMSFITS(self):
        """
        get rms fits according the calcode
        :return:
        """

        if self.calCode=="":
            return self.rmsRawCO12FITS


        if "12" in self.calCode and ("Raw" in self.calCode or "raw" in self.calCode):
            return self.rmsRawCO12FITS

        if "13" in self.calCode  and ("Raw" in self.calCode or "raw" in self.calCode):
            return self.rmsRawCO13FITS

        if "18" in self.calCode  and ("Raw" in self.calCode or "raw" in self.calCode):
            return self.rmsRawCO18FITS



        if "12" in self.calCode:
            return self.rmsCO12FITS

        if "13" in self.calCode:
            return self.rmsCO13FITS

        if "18" in self.calCode:
            return self.rmsCO18FITS



    def getSmoothRMS(self,outArmSmFITS):
        """

        :return:
        """

        #get all Outer Arm FITS, extract rms noise

        #smootFileList= self.getSmFITSFileList(calCode=self.codeLocalCO12)


        #testFile=smootFileList[0]
        saveName=self.addFixToFITSName( outArmSmFITS,"_rms" )

        saveName=self.rmsFITSPath+os.path.basename(saveName)

        doFITS.getRMSFITS(outArmSmFITS,saveName)

        return saveName
    def getRMSFITSByCodeAndsmf(self,calCode,smoothFactor):
        """

        :param calCode:
        :param smoothFactor:
        :return:
        """

        searchName="{}*SmFactor_{:.1f}*_rms.fits".format(calCode,smoothFactor)
        searchStr=self.rmsFITSPath+searchName

        searchResults=  glob.glob(searchStr)
        if len(searchResults) == 0 or len(searchResults) >1  :
            return None

        else:
            return searchResults[0]

    def isRaw(self,anyStr):
        """
        Check if "Raw" or "raw" in anyStr
        :param anyStr:
        :return:
        """
        return "Raw" in anyStr or "raw" in anyStr



    def selectByCode(self,calCode,selectList):



        searchStr = "12"

        if "12" in calCode:
            searchStr="12"

        if "13" in calCode:
            searchStr="13"

        if "18" in calCode:
            searchStr="18"

        isRawCode=self.isRaw(calCode)

        for eachElement in selectList:

            if  searchStr in eachElement and isRawCode==self.isRaw(eachElement):
                return eachElement


        return None






    def getOutArmCode(self,calCode):

        searchList=[self.codeOutCO12, self.codeOutCO13,self.codeOutCO18, self.codeRawOutCO12, self.codeRawOutCO13,self.codeRawOutCO18]
        return self.selectByCode(calCode,searchList)

    def getSmoothFactor(self,fitsFile):
        """

        :param smoothFactor:
        :return:
        """

        a,suffix=fitsFile.split("SmFactor_")

        try:
            firstFour = suffix[0:4]
            return np.float(firstFour)

        except:
            firstThree = suffix[0:3]
            return np.float( firstThree )

    def getNoiseFactor(self,fitsFile):
        """

        :param smoothFactor:
        :return:
        """

        a,suffix=fitsFile.split("NoiseAdd_")

        try:
            firstFour = suffix[0:3]
            return np.float(firstFour)

        except:
            firstThree = suffix[0:4]
            return np.float( firstThree )


    def checkRMSFITS(self,smFITS):
        """
        check the existance of the corresponding rms fits file, if not exist try to produce one
        :param smFITS:
        :return:
        """

        #step 1, search a file




        smFactor=self.getSmoothFactor(smFITS)

        outArmCode=self.getOutArmCode(self.calCode)


        rmsFITS=self.getRMSFITSByCodeAndsmf(outArmCode, smFactor)

        ########


        if rmsFITS==None:
            ###### produce the rms fits

            outArmCode=self.getOutArmCode(self.calCode)

            smFileOut=self.getSmFITSFileSingle(smFactor,calCode=outArmCode)


            if smFileOut is None:
                print "The corresponding Out arm smooth fits does not exist, could not produce rms fits, stop"
                return None

            else:
                #produce the rms fits with the out arm file
                #aaaaaaaa

                print "Out arm fits found, producing a rms fits accordingly."
                return  self.getSmoothRMS(smFileOut)
        else:
            print "The corresponding rms fits file produce with out arm fits found!"
            return rmsFITS



    def addAllNoise(self ):
        """
        find the smooth factor1, fits, add noise by absolut
        :return:
        """

        #first get all smooth fits


        allSmFiles=self.getSmFITSFileList()
        smFirstFile=allSmFiles[0]

        if "_SmFactor_1.0" not in smFirstFile:
            print "Wrong program, please check your data!"
            return


        for eachNFactor in self.noiseFactors:

            if eachNFactor==0.0:
                continue
            self.addNoiseByRMSFITS(smFirstFile,noiseFactor=eachNFactor)




    def addNoiseByRMSFITS(self,smFITS, noiseFactor=0.0, targetRMSFITS=None):

        
        """
        add noise spectral by spectral

        the smooth rms should be calculate with outArmFITS because the out arm has the largest about data points

        targetRMSFITS, should be the three data FITS files
        but the rawRMSFITS ,should be created with out arm fits cubes

        the noiseFactor is in Kelvin, which is the increase of noise to the smoothed FITS file

        :param smFITS:
        :param targetRmsFITS:
        :return:
        """

        saveSuffix = "{}{:.1f}".format(self.noiseTag,noiseFactor)
        saveName = self.addFixToFITSName(smFITS, saveSuffix, prefix=False)


        dataSM,headSM=doFITS.readFITS(smFITS)
        if "SmFactor_1.0" in smFITS and noiseFactor==0.0:
            #do not add noise to the raw fits file
            #raw data
            fits.writeto(saveName, dataSM, header=headSM, overwrite=True)
            doFITS.converto32bit(saveName, saveName)

            return saveName

        Nz,Ny,Nx=dataSM.shape

        #find corresponding Out Arm rms fits file, because out arm has the least signal range

        #outCode=self.getOutArmCode(self.calCode)
        #smFactor=self.getSmoothFactor(smFITS)

        if "SmFactor_1.0" in smFITS:

            rawRMSFITS=self.getRMSFITS() #use raw cofits


        else:
            rawRMSFITS=self.checkRMSFITS( smFITS  )


        if rawRMSFITS==None:
            print "Please check your data, the corresponding rmsFITS produce with out arm fits, does not exist"
            return

        rawRMSdata, rawRMShead =  doFITS.readFITS( rawRMSFITS )


        if targetRMSFITS==None:
            targetRMSFITS=self.getRMSFITS()

        targetRMSdata, targetRMShead = doFITS.readFITS( targetRMSFITS )

        targetRMSdata=noiseFactor+targetRMSdata # this step set the noise increse factor

        convolveRMS= np.sqrt( targetRMSdata**2-rawRMSdata**2 ) #any negative value data?


        print "the raw and the targed rms fits are"

        print rawRMSFITS,targetRMSFITS
        ##

        noiseData=self.produceNoiseFITS(Nz,inputRmsData=convolveRMS)

        print "Adding noise by {} K".format( noiseFactor )

        observedData = dataSM + noiseData


        fits.writeto( saveName , observedData,header= headSM,overwrite=True )
        doFITS.converto32bit(saveName,saveName)

        return saveName







    def getMinMeanMaxSize(self,TBList):


        minSizeList=[]
        meanSizeList=[]
        maxSizeList=[]

        for eachTB in TBList:

            size=self.getCloudSize(eachTB)

            minSizeList.append( np.min(size)   )
            meanSizeList.append(  np.mean(size)     )
            maxSizeList.append(  np.max(size)     )



        return minSizeList, meanSizeList, maxSizeList






    def drawCloudSizeChange(self,drawCode="max"):

        testCode=self.codeLocalCO12
        self.calCode = testCode
        TB,TBFiles=self.getAbsNoiseCleanTBList( absK=self.MWISPrmsCO12)
        #ceanTBFiless = self.getSmoothListFixNoise(getCleanTBFile=True)

        #print TBFiles

        minSizeList, meanSizeList, maxSizeList=self.getMinMeanMaxSize(TB)


        fig = plt.figure(figsize=(10, 8))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]


        axBeam = fig.add_subplot(1,1, 1)

        beamSize=self.smoothFactors*self.getBeamSize(self.calCode)

        #self.drawNoiseEffectSingle(axSens,1.5)
        if drawCode=="min":
            axBeam.plot(beamSize,minSizeList,'o-',color='b')
            axBeam.set_ylabel(r"Minimum cloud size (arcmin)")

        if drawCode == "max":
            axBeam.plot(beamSize, maxSizeList, 'o-', color='b')
            axBeam.set_ylabel(r"Maximum cloud size (arcmin)")

        if drawCode == "mean":
            axBeam.plot(beamSize, meanSizeList, 'o-', color='b')
            axBeam.set_ylabel(r"Mean cloud size (arcmin)")

        axBeam.set_xlabel(r"Beam Size (arcmin)")

        #self.drawColorBarNoiseFactor(axBeam)

        fig.tight_layout()
        plt.savefig("{}cloudNchange.pdf".format(drawCode) , bbox_inches='tight')
        plt.savefig("{}cloudNchange.png".format(drawCode), bbox_inches='tight', dpi=600)




    def readTBList(self,TBNameList):
        """

        :param TBNameList:
        :return:
        """
        #

        TBList=[]
        for eachName in TBNameList:
            TBList.append(Table.read(eachName) )
        return  TBList


    def drawCloudNumberChange(self):

        testCode=self.codeLocalCO12
        self.calCode = testCode
        #TB,TBFiles=self.getAbsNoiseCleanTBList( absK=self.MWISPrmsCO12)
        TBFiles = self.getSmoothListFixNoise(getCleanTBFile=True)
        TB= self.readTBList(TBFiles)

        #print TBFiles

        Nlist=self.getNList(TB)


        fig = plt.figure(figsize=(10, 8))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]


        axBeam = fig.add_subplot(1,1, 1)

        beamSize=self.smoothFactors*self.getBeamSize()

        #self.drawNoiseEffectSingle(axSens,1.5)

        #axBeam.plot(beamSize,Nlist,'o-',color='b')


        mwispF,cfaF,parameters=self.getFillingFactorAndDraw(  beamSize, Nlist,None, 1, drawFigure=False,saveTag="_Number")
        print mwispF,cfaF
        params,paramsErros=parameters
        axBeam.scatter(  beamSize  , Nlist , s=15,  color='red'   )

        fittingX = np.arange(0, np.max(beamSize), 0.01)
        # axFitting.plot( fittingX  ,  ffFunction(fittingX,params[0], params[1] , params[2]   ), color='blue'  )
        axBeam.plot(fittingX, ffFunction(fittingX, *params), color='blue', lw=1.5)


        axBeam.set_ylabel(r"Number of molecular clouds")
        axBeam.set_xlabel(r"Beam Size (arcmin)")

        #self.drawColorBarNoiseFactor(axBeam)
        axBeam.axvline(x=0, ls="--", color='black')

        fig.tight_layout()
        plt.savefig("cloudNchange.pdf", bbox_inches='tight')
        plt.savefig("cloudNchange.png", bbox_inches='tight', dpi=600)





    def getBeamSize(self ):
        if "12" in self.calCode:
            return self.rawBeamSizeCO12

        return self.rawBeamSizeCO13

    def getArmSubFTS(self):

        """
        Split fits into four arms
        :return:
        """
        #deakubg withCO 12
        CO12LocalAndOutFITS = "/home/qzyan/WORK/diskMWISP/fillingFactorData/q1Raw/Ua/mosaic_U.fits"
        CO12SgrAndscuFITS = "/home/qzyan/WORK/diskMWISP/fillingFactorData/q1Raw/Ub/mosaic_U.fits"

        CO13LocalAndOutFITS = "/home/qzyan/WORK/diskMWISP/fillingFactorData/q1Raw/La/mosaic_L.fits"
        CO13SgrAndscuFITS = "/home/qzyan/WORK/diskMWISP/fillingFactorData/q1Raw/Lb/mosaic_L.fits"

        CO18LocalAndOutFITS = "/home/qzyan/WORK/diskMWISP/fillingFactorData/q1Raw/L2a/mosaic_L2.fits"
        CO18SgrAndscuFITS = "/home/qzyan/WORK/diskMWISP/fillingFactorData/q1Raw/L2b/mosaic_L2.fits"



        doFITS.cropFITS(CO12LocalAndOutFITS, Vrange=self.outVrange,outFITS= self.outCO12FITS, overWrite=True )
        doFITS.cropFITS(CO12LocalAndOutFITS, Vrange=self.localVrange,outFITS=  self.localCO12FITS, overWrite=True )
        doFITS.cropFITS(CO12SgrAndscuFITS, Vrange=self.sgrVrange,outFITS= self.sgrCO12FITS, overWrite=True )
        doFITS.cropFITS(CO12SgrAndscuFITS, Vrange=self.scuVrange,outFITS= self.scuCO12FITS, overWrite=True )



        doFITS.cropFITS(CO13LocalAndOutFITS, Vrange=self.outVrange,outFITS= self.outCO13FITS, overWrite=True )
        doFITS.cropFITS(CO13LocalAndOutFITS, Vrange=self.localVrange,outFITS=  self.localCO13FITS, overWrite=True )
        doFITS.cropFITS(CO13SgrAndscuFITS, Vrange=self.sgrVrange,outFITS= self.sgrCO13FITS, overWrite=True )
        doFITS.cropFITS(CO13SgrAndscuFITS, Vrange=self.scuVrange,outFITS= self.scuCO13FITS, overWrite=True )



        doFITS.cropFITS(CO18LocalAndOutFITS, Vrange=self.outVrange,outFITS= self.outCO18FITS, overWrite=True )
        doFITS.cropFITS(CO18LocalAndOutFITS, Vrange=self.localVrange,outFITS=  self.localCO18FITS, overWrite=True )
        doFITS.cropFITS(CO18SgrAndscuFITS, Vrange=self.sgrVrange,outFITS= self.sgrCO18FITS, overWrite=True )
        doFITS.cropFITS(CO18SgrAndscuFITS, Vrange=self.scuVrange,outFITS= self.scuCO18FITS, overWrite=True )



    def checkCloudCubeSavePath(self ):
        """
        Examine the calCode of molecular clouds
        :param calCode:
        :return:
        """

        savePath= self.cloudCubePath+self.calCode


        if os.path.isdir( savePath ):
            pass
        else:
            os.mkdir(savePath)
        return savePath

    def getIndicesRaw(self, Z0, Y0, X0, values1D, choseID):

        cloudIndices = np.where(values1D == choseID)

        cX0 = X0[cloudIndices]
        cY0 = Y0[cloudIndices]
        cZ0 = Z0[cloudIndices]

        return tuple([cZ0, cY0, cX0])



    def addFFValues(self,fillingTB):
        """
        add filling factors to
        :return:
        """

        if type(fillingTB)==Table:
            ffTB=fillingTB
        else:

            ffTB = Table.read(fillingTB)
        ffTB = ffTB[ffTB[self.ffMWISPCol] >= 0]
        ffTB = ffTB[ffTB[self.ffMWISPCol] <= 1]
        ffTB = self.addMWISPFFerror(ffTB)
        ffTB = self.addCfaFFerror(ffTB)

        ffTB=ffTB[ ffTB[self.ffMWISPErrorCol]>0]


        return ffTB


    def fittingAngularSizeFF(self,fillingTB,showSizeRange=[0,40],saveTag="" ,useODR=False):
        """

        :return:
        """
        if saveTag=="":
            saveTag=self.calCode
        # get the distance splitting line between different

        drawCode = self.drawCodeSize

        #remove bad clouds

        ffTB= self.addFFValues(fillingTB)

        if 0: #do not use error control
            errorControl=ffTB[self.ffMWISPErrorCol]/ffTB[self.ffMWISPCol ]
            ffTB=ffTB[errorControl<0.3 ]


        print len(ffTB), "Number of total molecular clouds"

        pureVclipTB = self.pureVclip(ffTB)

        pureLBclipTB = self.pureLBclip(ffTB)

        noClipTB = self.noClipClouds(ffTB)
        bothClipTB = self.bothVAndLBClip(ffTB)

        edgeClodus = self.getEdgeClouds(ffTB)

        print "Number of edge clouds", len(edgeClodus)

        print  len(ffTB) - len(noClipTB) - len(pureLBclipTB) - len(pureVclipTB) - len(bothClipTB), "??????"

        fig = plt.figure(figsize=(15, 6))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})

        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        
        minVel = np.min(ffTB["v_cen"])
        maxVel = np.max(ffTB["v_cen"])

        # axFF.scatter(   ffTB["area_exact"]/3600., ffTB[self.ffMWISPCol]  , color="blue",  s=3)
        # ffMWISPCol   ffCfACol

        elinewidth = 0.6
        markerSize = 2.1

        # select those molecular clouds, with filling factors less than 0.3
        #noClipTB = noClipTB[noClipTB[self.ffMWISPCol] < 0.3]

        #use all clouds
        ################################################# fitting data
        useTB= ffTB
        
        drawX = self.getCloudSize(useTB)
        drawY = useTB[self.ffMWISPCol]
        yError= useTB[self.ffMWISPErrorCol]
        Vs = useTB["v_cen"]

        #radius= drawX/2
        #perimeterPixN=2*np.pi*radius/0.5

        errorArea =   np.sqrt( useTB["area_exact"]/0.25 )*0.25

        errorSize= errorArea/np.sqrt(useTB["area_exact"])/np.sqrt(np.pi)

        print np.min(  useTB["area_exact"]),"miniumum area"
        #####################################################
        cmap = plt.cm.jet
        # norm = matplotlib.colors.BoundaryNorm(np.arange(0,30,0.1), cmap.N)
        #normV = mpl.colors.Normalize(vmin=np.min(Vs), vmax=np.max(Vs))
        normV = mpl.colors.Normalize(vmin=np.min(Vs), vmax=np.max(Vs) )

        m = plt.cm.ScalarMappable(norm=normV, cmap=cmap)
        v_color = np.array([(m.to_rgba(v)) for v in Vs])


        #draw fffactor with velocity colors
        #######################################################################################################
        axFFvel =  fig.add_subplot(1, 2, 1 )

        #sc = self.drawErrorBar(axFFvel, noClipTB, drawCode=drawCode, markerSize=markerSize, color='gray', markerType=".",   elinewidth=elinewidth, label="Complete in PPV space", showYError=False)
        # self.drawErrorBar(axFF,pureVclipTB, drawCode =drawCode, markerSize=markerSize+0.8,color='b',markerType="D",elinewidth=elinewidth,label="Incomplete in v space" ,showYError=False)
        # self.drawErrorBar(axFF,pureLBclipTB, drawCode =drawCode, markerSize=markerSize+0.8,color='r',markerType="^",elinewidth=elinewidth,label="Incomplete in l-b space" ,showYError=False)
        sc = axFFvel.scatter(drawX, drawY, c=Vs, cmap=cmap, norm=normV, s=3, facecolors='none', lw=0.5, marker="o", label="")

        axFFvel.set_xlim( showSizeRange )

        if 1:
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(axFFvel)
            cax1 = divider.append_axes("right", size="3%", pad=0.05)
            cb = fig.colorbar(sc, cax=cax1)
            cax1.set_ylabel(r"$V_{\rm LSR}$ ($\rm km\ s^{-1}$)")
        ###########################################################################################################
        axFF = fig.add_subplot(1, 2, 2,sharex=axFFvel,sharey=axFFvel)

            # c =Vs, cmap=cmap,norm=normV,
            # ax.errorbar(drawX,   drawY , yerr= yerr ,    markersize=markerSize , linestyle='none', c= v_color ,     marker= markerType , capsize=0,  elinewidth=elinewidth , lw=1, label= label )
            # ax.errorbar(drawX,   drawY , yerr= yerr ,    markersize=markerSize , linestyle='none', c= Vs ,  cmap=cmap,norm=normV,   marker= markerType , capsize=0,  elinewidth=elinewidth , lw=1, label= label )
            #sc = ax.scatter(drawX, drawY, c='gray', cmap=cmap, norm=normV, s=3, facecolors='none', lw=0.5, marker="o",  label="")
        #axFF.scatter(drawX, drawY, c="gray", cmap=cmap,   s=3, facecolors='none', lw=0.5, marker="o",label="")
        #axFF.errorbar(drawX,   drawY , xerr=errorSize, yerr= yError ,    markersize=1 , fmt='o',  color='gray' , capsize=0.1,  elinewidth=0.5 , lw=1 ,alpha=0.8  ,label="" )
        axFF.errorbar(drawX,   drawY ,  yerr= yError ,    markersize=1 , fmt='o',  color='gray' , capsize=0.1,  elinewidth=0.5 , lw=1 ,alpha=0.8  ,label="" )

        angularSize1 = 10.5
        #axFF.axvline(x=angularSize1, ls="--", color='black')

        # find the mean Veloicty, of angularSize1

        sub1 = self.clipByAngularSize(noClipTB, angularSize1, 13)
        print np.mean(sub1["v_cen"]), np.std(sub1["v_cen"], ddof=1)

        angularSize2 = 9.0
        #axFF.axvline(x=angularSize2, ls="--", color='black')
        sub2 = self.clipByAngularSize(noClipTB, angularSize2, angularSize1)
        print np.mean(sub2["v_cen"]), np.std(sub2["v_cen"], ddof=1)

        # search  velocity layers, based on noClipTB

        #axFF.legend(loc=4)

        ####################
        #show the function of overall filling factor and angular size relationship



        at = AnchoredText(saveTag, loc=1, frameon=False)
        axFF.add_artist(at)

        axFF.set_ylim([-0.1, 1.1])

        axFFvel.set_ylabel("Filling factor")
        axFF.set_xlabel("Angular size (arcmin)")
        axFFvel.set_xlabel("Angular size (arcmin)")


        #fitting the parameters
        ###################################################
        #print errorSize

        para,paraError=self.fittingFFAndSize(testffAndSizeFunc, drawX,drawY, yError,sizeError=errorSize,  useODR=useODR)
        #para,paraError=self.fittingFFAndSize(ffFunction , drawX,drawY,  yError)

        #para,paraError=self.ffOdr( drawX,drawY,yError, xError=None )
        #para,paraError=self.ffOdr( drawX,drawY,yError, xError=errorSize )

        #para,paraError=self.fittingFFSizeRelationOdr(  drawX, drawY, None, yError)
        #para,paraError=self.fittingFFSizeRelationOdr(  drawX, drawY, errorSize, yError)


        #using MCMC, looks like very little difference with
        #from MCMCfunctionFitting import fittingFFMCMC
        #doFFmcmc=fittingFFMCMC()
        #para,paraError=doFFmcmc.getParameters(  drawX,drawY,  yError)





        print "The fitting result is "

        print para
        print paraError
        if len(para)==2:
            a,b=para

            x=np.arange(-0.1,showSizeRange[1],0.01)
            formulaTex=r"$y={:.3f}-\frac{{{:.3f}}}{{x+{:.3f}}}$".format(a,b,b/a)
            axFF.plot(x, testffAndSizeFunc(x,a,b),  "b-",lw=1,label= formulaTex )

            #add formula in the figure
            axFF.legend(loc=4)

            zeroPoint = 0 #0.564  # 0.977

            #c = b / a - zeroPoint


            axFF.axvline(x=zeroPoint,ls="--",color='black',lw=0.8)

        if  len(para)==3:
            a,b,c=para
            x=np.arange(-0.1,showSizeRange[1],0.01)
            formulaTex=r"$y={:.3f}-\frac{{{:.3f}}}{{x+{:.3f}}}$".format(a,b,c)
            axFF.plot(x, testffAndSizeFunc(x,a,b,c),  "b-",lw=1,label= formulaTex )

            #add formula in the figure
            axFF.legend(loc=4)
            print "ZeropPoint???????????????????????/",b/a-c

        if  len(para)==1:
            a =para[0]
            x=np.arange(-0.1,showSizeRange[1],0.01)
            formulaTex=r"$y= \frac{{x}}{{x+{:.3f}}}$".format(a )
            axFF.plot(x, testffAndSizeFunc(x,a  ),  "b-",lw=1,label= formulaTex )

            #add formula in the figure
            axFF.legend(loc=4)



            axFF.axvline(x=0,ls="--",color='black',lw=0.8)

        axFF.axhline(y=0,ls="--",color='green',lw=0.8)
        plt.subplots_adjust(wspace=0.3)
        axFF.set_xlim(showSizeRange )
        plt.savefig(self.paperFigurePath+"{}FittingAngularSizeRelation.png".format(saveTag), bbox_inches='tight', dpi=300)



    def testDistanceLayers(self, ):
        """

        :return:
        """

        #get the distance splitting line between different

        drawCode=self.drawCodeSize
        #draw the


        fillingTB =  "edgeInfo_fillingFactor_fluxTB_rawLocalCO12rawLocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"


        ffTB= Table.read( fillingTB )
        ffTB=ffTB[ffTB[self.ffMWISPCol]>0]

        ffTB=self.addMWISPFFerror(ffTB)
        ffTB=self.addCfaFFerror(ffTB)





        print len(ffTB),"Number of total molecular clouds"

        pureVclipTB=self.pureVclip(ffTB)

        pureLBclipTB = self.pureLBclip( ffTB)

        noClipTB= self.noClipClouds(ffTB)
        bothClipTB = self.bothVAndLBClip(ffTB)

        edgeClodus= self.getEdgeClouds(ffTB)

        print "Number of edge clouds", len(edgeClodus)


        print  len(ffTB)- len(noClipTB ) -  len(pureLBclipTB )-  len(pureVclipTB ) - len(bothClipTB) ,"??????"


        fig = plt.figure(figsize=(13, 6))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})


        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        axFF = fig.add_subplot(1,2, 1)


        minVel = np.min(ffTB["v_cen"])
        maxVel = np.max(ffTB["v_cen"])

        #axFF.scatter(   ffTB["area_exact"]/3600., ffTB[self.ffMWISPCol]  , color="blue",  s=3)
        # ffMWISPCol   ffCfACol

        elinewidth = 0.6
        markerSize=2.1


        #select those molecular clouds, with filling factors less than 0.3
        noClipTB=noClipTB[noClipTB[self.ffMWISPCol]<0.3]


        if drawCode==self.drawCodeSize: #angular size


            sc=self.drawErrorBar(axFF,noClipTB, drawCode =drawCode, markerSize=markerSize,color='gray',markerType=".",elinewidth=elinewidth,label="Complete in PPV space" ,showYError=False)
            #self.drawErrorBar(axFF,pureVclipTB, drawCode =drawCode, markerSize=markerSize+0.8,color='b',markerType="D",elinewidth=elinewidth,label="Incomplete in v space" ,showYError=False)
            #self.drawErrorBar(axFF,pureLBclipTB, drawCode =drawCode, markerSize=markerSize+0.8,color='r',markerType="^",elinewidth=elinewidth,label="Incomplete in l-b space" ,showYError=False)

            axFF.set_xlim([0,15 ])
            axFF.set_xlabel("Angular size (arcmin)")

            if 1:
                from mpl_toolkits.axes_grid1 import make_axes_locatable
                divider = make_axes_locatable(axFF)
                cax1 = divider.append_axes("right", size="3%", pad=0.05)
                cb = fig.colorbar(sc, cax=cax1)
                cax1.set_ylabel(r"$V_{\rm LSR}$ ($\rm km\ s^{-1}$)")

        #draw AngularSizeLines

        angularSize1=10.5
        axFF.axvline(x= angularSize1 , ls="--", color='black')

        #find the mean Veloicty, of angularSize1

        sub1=self.clipByAngularSize(noClipTB,angularSize1,13 )
        print np.mean(sub1["v_cen"]), np.std(sub1["v_cen"],ddof=1)



        angularSize2=9.0
        axFF.axvline(x= angularSize2 , ls="--", color='black')
        sub2=self.clipByAngularSize(noClipTB,angularSize2,angularSize1 )
        print np.mean(sub2["v_cen"]), np.std(sub2["v_cen"],ddof=1)



        #search  velocity layers, based on noClipTB

        axFF.legend(loc=4)

        #draw
        #### add Region to the figure


        at = AnchoredText(self.calCode , loc=1, frameon=False)
        axFF.add_artist(at)


        axFF.set_ylim([-0.2, 1.2  ])


        axFF.set_ylabel("Filling factor")

        ################################## axvel
        axVel = fig.add_subplot(1,2, 2)

        angularSizeDraw=self.getCloudSize(noClipTB)

        axVel.scatter(angularSizeDraw,  noClipTB["v_cen"], s=5)
        #axVel.invert_xaxis()
        plt.subplots_adjust(  wspace=0.4)
        axVel.set_ylabel("VLSR (km/s)")
        axVel.set_xlabel("Angular size (arcmin)")


        plt.savefig( "{}distancelayers_{}.png".format(self.calCode, drawCode ), bbox_inches='tight', dpi=600)

    def clipByAngularSize(self,TB,minAng,maxAng):

        angularSizeAll=self.getCloudSize(TB)

        selectionCriteria=np.logical_and( angularSizeAll>=minAng, angularSizeAll<=maxAng  )
        return TB[selectionCriteria]



    def getCloudCubes(self,  labelsFITS,cloudTBFile , calCode=None, writeFITS=False, saveTag="" ):

        """

        #output all data cubes for each cloud

        :return:
        """

        #################

        if calCode!=None:

            self.calCode=calCode
        savePath =  self.checkCloudCubeSavePath() #"./cloudSubCubes/"

        print "Processing fits Label ", labelsFITS
        print "Processing TB file", cloudTBFile

        rawCOFITS=self.getRawCOFITS()
        rmsFTS=self.getRMSFITS()

        cloudTB = Table.read(cloudTBFile)
        dataCluster, headCluster = myFITS.readFITS(labelsFITS)
        dataCO, headCO = myFITS.readFITS( rawCOFITS)
        rmsData,rmsHead=myFITS.readFITS( rmsFTS )

        _,Vaxis =  doFITS.getSpectraByIndex(dataCO, headCO, 0, 0)
        #consecutiveMask=self.getConsecutiveMask(dataCO,rmsData,cutoff=3, conNumber=3)



        minV = np.nanmin(dataCluster[0])
        wcsCloud = WCS(headCluster)
        clusterIndex1D = np.where(dataCluster > minV)
        clusterValue1D = dataCluster[clusterIndex1D]
        Z0, Y0, X0 = clusterIndex1D

        fitsZero = np.zeros_like(dataCluster,dtype=np.float)


        #add a function that check if the cloud touches L edge,B edges, or zEdges

        Nz, Ny, Nx = dataCO.shape
        zeroProjection = np.zeros((Ny, Nx))  # one zero channel, used to get the projection area and
        zeroProjectionSpectra = np.zeros((Ny, Nx))  # one zero channel, used to get the projection area and

        ###### if this cloud touches the edge
        cloudTB[self.touchLedgeCol]  = cloudTB["peak"]*0
        cloudTB[self.touchLedgeCol].unit=""

        
        cloudTB[self.touchBedgeCol]  =cloudTB["peak"]*0
        cloudTB[self.touchBedgeCol].unit=""

        
        cloudTB[self.touchVedgeCol]  =  cloudTB["peak"]*0
        cloudTB[self.touchVedgeCol].unit=""

        cloudTB[self.peakChannelCol] = cloudTB["peak"]*0
        cloudTB[self.peakChannelCol].unit=""

        cloudTB[self.conseCol] = cloudTB["peak"]*0 #number of spectra that have
        cloudTB[self.conseCol].unit=""
        #cloudTB["peakB"].unit=""
        #cloudTB["peakL"].unit=""
        #cloudTB["peakV"].unit=""
        #cloudTB["peak"].unit="K"

        cloudTB[self.eqLwCol] = cloudTB["peak"]*0 #number of spectra that have
        cloudTB[self.eqLwCol].unit="km/s"

        cloudTB[self.vrmsAvgSpectraCol] = cloudTB["peak"]*0 #number of spectra that have
        cloudTB[self.vrmsAvgSpectraCol].unit="km/s"

        cloudTB[self.vmeanAvgSpectraCol] = cloudTB["peak"]*0 #number of spectra that have
        cloudTB[self.vmeanAvgSpectraCol].unit="km/s"

        cloudTB.sort("allChannel")
        #check, the spectral of the peak position, to get the number channel, that the spectra toward the peak position must have the
        # three consecutive channels

        #plot figures
        fig = plt.figure(figsize=(10, 6  ))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 16, 'serif': ['Helvetica']})


        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]


        axSpec = fig.add_subplot(1,1, 1)

        drawFigure=False

        widgets = ['Extracting cloud cubes: ', Percentage(), ' ', Bar(marker='0', left='[', right=']'), ' ', ETA(), ' ',
                   FileTransferSpeed()]  # see docs for other options
        pbar = ProgressBar(widgets=widgets, maxval=len(cloudTB))
        pbar.start()
        i=0

        #do not use over all spectra, use the peak spectra

        for eachC in cloudTB:
            plt.cla()
            i=i+1

            cloudID = eachC["_idx"]
            saveName = "{}cloud{}cube.fits".format(self.calCode,cloudID)

            cloudIndex = self.getIndicesRaw(Z0, Y0, X0, clusterValue1D, cloudID)
            cloudZ0, cloudY0, cloudX0 = cloudIndex
            projectIndex = tuple([cloudY0, cloudX0])
            zeroProjection[projectIndex] = 1

            #getAverageSpectra, and

            #fitsZero[testIndices] = dataCO[testIndices]
            fitsZero[cloudIndex] =  dataCO[cloudIndex]
            if 1: #calculate velocity dispersion, two versions
                # cropThe cloudRange
                minY = np.min(cloudY0)
                maxY = np.max(cloudY0)

                ###########
                minX = np.min(cloudX0)
                maxX = np.max(cloudX0)

                ###########
                minZ = np.min(cloudZ0)
                maxZ = np.max(cloudZ0)
                ######################################################################
                #cloudCropSpectra = fitsZero[:, minY:maxY + 1, minX:maxX + 1]

                cloudCropCube = fitsZero[minZ:maxZ + 1, minY:maxY + 1, minX:maxX + 1]
                cloudVaxis= Vaxis[ minZ:maxZ+1 ]
                
                averageSpectraCrop = np.nansum(cloudCropCube, axis=(1, 2))

                #intCloud = np.nansum(cloudCropCube, axis=0)

                # count the number spectra

                totalSpectral = len(zeroProjection[zeroProjection > 0])

                meanSpectral = averageSpectraCrop / 1. / totalSpectral

                spectraPeak = np.max(meanSpectral)

                area = self.getVelResolution() * np.sum(meanSpectral)

                eqLineWidth = area / spectraPeak
                eachC[self.eqLwCol] = eqLineWidth  # np.max( [totalPeakN1, totalPeakN2, totalPeakN3] )


                vMeanAverageSp,vrmsAverageSp=doFITS.weighted_avg_and_std(cloudVaxis, meanSpectral   )
                eachC[self.vrmsAvgSpectraCol] = vrmsAverageSp  # np.max( [totalPeakN1, totalPeakN2, totalPeakN3] )

                eachC[self.vmeanAvgSpectraCol] = vMeanAverageSp  # np.max( [totalPeakN1, totalPeakN2, totalPeakN3] )

            #fitsZero[cloudIndex] =consecutiveMask[cloudIndex]  #too slow, need a better way
            #consecutiveCheck= consecutiveMask[cloudIndex]

            #goodSpectraIndex = self.getIndicesRaw(Z0, Y0, X0, consecutiveCheck, 1)
            #spectraZ0, spectraY0, spectraX0 = goodSpectraIndex
            #spectraPosition = tuple([spectraY0, spectraX0])

            #zeroProjectionSpectra[spectraPosition]=1


            #spectraCloud=np.nanmax(fitsZero,axis=0)

            #eachC[self.conseCol]= np.sum( zeroProjectionSpectra  ) #np.nanmax( consecutiveCheck) #how many spectra has
            #zeroProjectionSpectra[spectraPosition]=0



            #spectraN= eachC["area_exact"]/0.25 #np.sum(zeroProjection)

            #spectraIndex=np.where(zeroProjection==1)

            #spIndexY,spIndexX=spectraIndex

            #spSum=dataCO[:,0,0]*0

            #for j in range(len(spIndexX)):

                #spB=spIndexY[j]
                #spL= spIndexX[j]
                #spSum=spSum+dataCO[:,spB,spL]


            #calculate averageSpectra

            #spectraSum=dataCO*zeroProjection
            #print "Equal?",spectraN,np.sum(zeroProjection)
            #averageSpectra =  spSum /spectraN


            zeroProjection[projectIndex] = 0 #untag zeroProjection

            minZ = np.min(cloudZ0)
            maxZ = np.max(cloudZ0)

            minY = np.min(cloudY0)
            maxY = np.max(cloudY0)

            minX = np.min(cloudX0)
            maxX = np.max(cloudX0)

            ####
            if minZ==0 or maxZ==Nz-1:
                eachC[self.touchVedgeCol]=1


            if minY==0 or maxY==Ny-1:
                eachC[self.touchBedgeCol]=1

            if minX==0 or maxX==Nx-1:
                eachC[self.touchLedgeCol]=1

            #record the peakSpectra
            peakB= int( eachC["peakB"] ) #index

            peakL= int(  eachC["peakL"] )

            peakV=int(eachC["peakV"])

            #cloudPeakSpectra,vel=doFITS.getSpectraByIndex(fitsZero,headCluster,peakL,peakB)


            #use average spectral to detect bad channels, the method is to get an average spectra, find the good channel number
            #should be 3*3
            if 0:#average spectra of raw CO without soMask
                #averageSize=3

                #lowBcut=max([0,peakB-averageSize/2])
                #upBcut=min([Ny,peakB+averageSize])
                #lowLcut=max([0,peakL-averageSize/2])
                upLcut=min([Nx,peakL+averageSize])


                cropData=dataCO[:,lowBcut:upBcut,lowLcut:upLcut]

                averageSpectra=np.nanmean(cropData,axis=(1,2))
                #calculate rms
                cutSigma= 3
                negativeData=averageSpectra[averageSpectra<0]
            #save figures
            if drawFigure:
                axSpec.step(range(Nz),averageSpectra,'b-',lw=0.5,where='mid')



            #if len(negativeData)<10:
                #eachC[self.peakChannelCol] = Nz-10 # np.max( [totalPeakN1, totalPeakN2, totalPeakN3] )
            if 1:
                pass

            else:
                rmsSpectra=np.std(negativeData,ddof=1)/np.sqrt(1-2./np.pi )
                cloudPeakSpectra=averageSpectra>= cutSigma*rmsSpectra
                cloudPeakSpectra=cloudPeakSpectra*1

                allStd=np.nanstd(averageSpectra,ddof=1)
                if drawFigure:
                    axSpec.axvline(x=peakV, ls="--", color='black',lw=0.5,label="Spectral sigma:{:.2f}, overall std{:.2f}".format(rmsSpectra,allStd) )

                if cloudPeakSpectra[peakV] < 1:
                    eachC[self.peakChannelCol] = 0  #  #peak less then 3 sigma, reject
                    if drawFigure:

                        axSpec.axhline(y=cutSigma * rmsSpectra, ls="--", color='red', lw=0.8,  label="Channel number above {} sigma: {}".format(cutSigma, 0))


                else:
                    totalN=1

                    moveIndex=1 #move foreward
                    while 1:
                        if peakV+moveIndex>=Nz:
                            break
                        if  cloudPeakSpectra[peakV+moveIndex ]<1:
                            break

                        totalN=totalN+1
                        moveIndex=moveIndex+1

                    moveIndex=1 #move  backword
                    while 1:
                        if peakV-moveIndex<0 :# or cloudPeakSpectra[peakV-moveIndex ]<1 :
                            break

                        if  cloudPeakSpectra[peakV-moveIndex ]<1 :
                            break

                        totalN=totalN+1
                        moveIndex=moveIndex+1
                    if drawFigure:
                        axSpec.axhline(y=cutSigma* rmsSpectra, ls="--", color='red', lw=0.8,  label="Channel number above {} sigma: {}".format(cutSigma, totalN))

                    #find the consecutive number channels around the cloud peak
                    #totalPeakN1=   np.sum(  cloudPeakSpectra[peakV-2:peakV+1] )
                    #totalPeakN2= np.sum( cloudPeakSpectra[peakV-1:peakV+2] )
                    #totalPeakN3= np.sum(  cloudPeakSpectra[peakV:peakV+3] )

                    #this should be three
                    eachC[self.peakChannelCol] =  totalN  #np.max( [totalPeakN1, totalPeakN2, totalPeakN3] )
            if drawFigure:

                axSpec.legend()
                plt.savefig(self.badChannelSpecPath+"{}Cloud_{}_Spec.png".format(self.calCode, cloudID), bbox_inches='tight', dpi=100)





            if writeFITS:
                cropWCS = wcsCloud[minZ:maxZ + 1, minY:maxY + 1, minX:maxX + 1]
                saveFullName = os.path.join(savePath, saveName)
                #fitsZero[cloudIndex] =  dataCO[cloudIndex]

                cropData = fitsZero[minZ:maxZ + 1, minY:maxY + 1, minX:maxX + 1]

                fits.writeto(saveFullName , cropData, header=cropWCS.to_header(), overwrite=True)

            fitsZero[cloudIndex] = 0

            pbar.update(i)
        pbar.finish()

        cloudTB["sum"].unit="K"
        #del cloudTB["area_ellipse"]
        #del cloudTB["flux"]
        #del cloudTB["major_sigma"]
        #del cloudTB["minor_sigma"]
        #del cloudTB["radius"]

        cloudTB["v_cen"].unit="km/s"
        cloudTB["x_cen"].unit="deg"
        cloudTB["y_cen"].unit="deg"

        
        if saveTag=="":

            cloudTB.write("edgeInfo_"+os.path.basename(cloudTBFile),overwrite=True)
            cloudTBRemoveFakeCloud=self.removeFakeClouds(cloudTB)

            cloudTBRemoveFakeCloud.write("pureEdgeInfo_"+os.path.basename(cloudTBFile),overwrite=True)

            

        else:
            cloudTB.write(saveTag+os.path.basename(cloudTBFile),overwrite=True)


    def getSelectCriteria(self,TB,colName,lRange):
        """

        :param TB:
        :param colName:
        :param range:
        :return:
        """

        minL = min(lRange)
        maxL = max(lRange)

        selectLPart1 = TB[colName] >= minL
        selectLPart2 = TB[colName] <= maxL

        selectL = np.logical_and(selectLPart1, selectLPart2)
        return selectL

    def removeByLBVrange(self,TB,lRange=None,bRange=None,vRange=None):
        """

        :param TB:
        :param lRange:
        :param bRange:
        :param vRange:
        :return:
        """



        badCloudSelect=TB[self.VCol] ==TB[self.VCol]

        if lRange!=None:

            selectL= self.getSelectCriteria(TB,self.LCol,lRange)
            badCloudSelect=np.logical_and(badCloudSelect,selectL)

        if bRange!=None:

            selectB= self.getSelectCriteria( TB, self.BCol, bRange)
            badCloudSelect=np.logical_and( badCloudSelect, selectB )


        if vRange!=None:

            selectV= self.getSelectCriteria( TB, self.VCol, vRange)
            badCloudSelect=np.logical_and( badCloudSelect, selectV )




        return TB[~badCloudSelect] #only return good clouds





    def removeFakeClouds(self,TB):
        """
        need to know the calcode, otherwise, do not know which line
        :param TB:
        :return:
        """

        newTB=TB.copy()

        if self.calCode in self.rawLocalCodeList: #Local clouds are all good
            return newTB


        if self.calCode==self.codeRawOutCO12: #Outer Arm, CO12
            vRange1=[-42,-41]

            lRange,bRange=doFITS.box(45.1168245, -2.9644784, 31924.800,14860.800, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange1,lRange=lRange,bRange=bRange)

            lRange,bRange=doFITS.box(37.0309402, -2.4327191, 11729.167,13500.000, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange1,lRange=lRange,bRange=bRange)

            lRange,bRange=doFITS.box(42.7673245, 3.4924883, 17070.000,11250.000, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange1,lRange=lRange,bRange=bRange)

            lRange,bRange=doFITS.box(39.6126649, 1.5151520, 2838.260,1314.475, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange1,lRange=lRange,bRange=bRange)

            lRange,bRange=doFITS.box(37.4773628, 0.6109141, 1995.434,760.497, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange1,lRange=lRange,bRange=bRange)

            lRange,bRange=doFITS.box(37.3996547, 1.6397105, 5310.000,3720.000, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange1,lRange=lRange,bRange=bRange)

            lRange,bRange=doFITS.box(34.7809529, 3.7184850, 26948.171,10151.389, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange1,lRange=lRange,bRange=bRange)

            lRange, bRange =doFITS.box(26.6407364, 0.4167354, 2879.051,1837.384, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange1,lRange=lRange,bRange=bRange)

            lRange,bRange=doFITS.box(27.4384634, 1.4937673, 1678.241,737.847, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange1,lRange=lRange,bRange=bRange)

            lRange,bRange=doFITS.box(26.0317835, 1.7200975, 1589.821,288.803, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange1,lRange=lRange,bRange=bRange)


            print "????????????????????????"

            vRange2=[-25,-23.5]

            lRange,bRange=doFITS.box(47.2802382, 3.2289989, 17050.000,14225.000, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange2,lRange=lRange,bRange=bRange)

            lRange,bRange=doFITS.box(44.8463726, -3.4107014, 19592.233,17189.190, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange2,lRange=lRange,bRange=bRange)

            lRange,bRange=doFITS.box(31.9797976, -2.7071768, 14228.500,18015.167, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange2,lRange=lRange,bRange=bRange)

            lRange, bRange =doFITS.box(32.4924728, 1.4271889, 2220.374,5927.694, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange2,lRange=lRange,bRange=bRange)

            lRange,bRange=doFITS.box(30.6750868, 1.1630478, 2501.688,4300.090, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange2,lRange=lRange,bRange=bRange)

            lRange,bRange=doFITS.box(30.9622206, 0.2806380, 1896.286,651.433, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange2,lRange=lRange,bRange=bRange)

            return newTB


        if self.calCode==self.codeRawOutCO13: #Outer Arm, CO13

            newTB=self.removeByLBVrange(newTB, vRange=[-75,-70] )
            newTB=self.removeByLBVrange(newTB, bRange=[3, 5] )

            return newTB



        if self.calCode==self.codeRawOutCO18: #Outer Arm, CO18

            newTB=self.removeByLBVrange(newTB, vRange=[-27,-21] )
            newTB=self.removeByLBVrange(newTB, vRange=[-71,-67] )



            return newTB




        if self.calCode==self.codeRawScuCO12: #ScuTum arm CO12

            newTB=self.removeByLBVrange(newTB, vRange=[115,117], lRange=[34.5, 50] )
            newTB=self.removeByLBVrange(newTB, vRange=[130, 135 ], lRange=[26.25,50] )

            return newTB

        if self.calCode==self.codeRawScuCO13: #ScuTum arm CO13



            vRange1=[95,98]


            lRange, bRange = doFITS.box(43.8741486,0.0333496,44149.238,37929.792,0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange1,lRange=lRange,bRange=bRange)

            lRange, bRange = doFITS.box(37.2196234,-2.6091924,6882.940,19481.889,0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange1,lRange=lRange,bRange=bRange)

            lRange, bRange = doFITS.box(31.1818432,-3.2707571,39225.600,14385.600,0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange1,lRange=lRange,bRange=bRange)

            lRange, bRange = doFITS.box(33.4215766,3.3910429,33177.600,13115.520,0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange1,lRange=lRange,bRange=bRange)

            lRange, bRange = doFITS.box(36.7022678,1.2874327,4575.803,4873.495,0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange1,lRange=lRange,bRange=bRange)

            lRange, bRange = doFITS.box(27.4948766,1.5292040,2950.000,2725.000,0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange1,lRange=lRange,bRange=bRange)



            vRange2=[74,76]
            lRange, bRange = doFITS.box(38.5967579,4.7255613,2930.357,3215.021,0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange2,lRange=lRange,bRange=bRange)



            return newTB
        if self.calCode==self.codeRawScuCO18: #ScuTum arm CO13

            vRange1=[112, 114]






            lRange, bRange = doFITS.box(30.6786667,-3.1768667,24192.000,14385.600,0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange1,lRange=lRange,bRange=bRange)

            lRange, bRange = doFITS.box(45.0777176,-0.0616218,37138.789,37858.429,0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange1,lRange=lRange,bRange=bRange)


            vRange2=[92, 98]

            lRange, bRange = doFITS.box(43.9179707, -0.1222090, 43291.447,37935.241, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange2,lRange=lRange,bRange=bRange)

            lRange, bRange = doFITS.box(35.1288007, 2.0227583, 19124.980,10812.488, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange2,lRange=lRange,bRange=bRange)

            lRange, bRange = doFITS.box(37.2892114, -2.6090439, 9151.631,18906.537, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange2,lRange=lRange,bRange=bRange)

            lRange, bRange = doFITS.box(37.5133908, 0.1635380, 712.240,2088.269, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange2,lRange=lRange,bRange=bRange)

            lRange, bRange = doFITS.box(33.1570346, -3.2696701, 26000.052,14000.028, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange2,lRange=lRange,bRange=bRange)

            lRange, bRange = doFITS.box(36.8156083, 0.2946295, 831.888,759.550, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange2,lRange=lRange,bRange=bRange)



            vRange3=[93.5, 95]

            lRange, bRange = doFITS.box(28.1729102, 0.4796709, 1624.870,2450.886, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange3,lRange=lRange,bRange=bRange)

            lRange, bRange = doFITS.box(27.8044150, 0.4677465, 978.973,1656.463, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange3,lRange=lRange,bRange=bRange)



            return  newTB

        if self.calCode==self.codeRawSgrCO12: #Sagitarrius arm CO12
            return  newTB

        if self.calCode==self.codeRawSgrCO13: #Sagitarrius arm CO13
            newTB=self.removeByLBVrange(newTB, vRange=[56,58],lRange=[45,50],bRange= [2,5] )

            return newTB

        if self.calCode == self.codeRawSgrCO18:  #Sagitarrius arm CO13

            vRange1=[60,62]

            lRange, bRange = doFITS.box(44.6240614, 1.9979111, 23675.968,11345.378, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange1,lRange=lRange,bRange=bRange)

            lRange, bRange = doFITS.box(39.6437782, 0.8537092, 3094.234,6806.678, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange1,lRange=lRange,bRange=bRange)

            lRange, bRange = doFITS.box(38.0645402, 0.5305706, 1319.444,2968.750, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange1,lRange=lRange,bRange=bRange)

            lRange, bRange = doFITS.box(38.8616970, 0.0018808, 1151.212,1976.828, 0)
            newTB=self.removeByLBVrange(newTB, vRange=vRange1,lRange=lRange,bRange=bRange)


            vRange2=[61,63]

            lRange, bRange = doFITS.box(40.9764480, 4.9925406, 3298.611,1215.278, 0)

            newTB=self.removeByLBVrange(newTB, vRange=vRange2,lRange=lRange,bRange=bRange)




            return newTB
        
        
        
        
        

        return newTB


    def getTotalFlux(self,fitsFile):
        """
        Use DBSCAN to maked to
        :param fitsFile:
        :return:
        """




    def smoothFITSbySMFactor(self,rawCOFITS=None, rawBeamSize = 52./60  ):

        """
        Default
        :param rawCOFITS:
        :param rawBeamSize:
        :return:
        """
        #calCode=rawCOFITS[0:-5]
        #self.calCode=calCode

        if rawCOFITS is None:
            rawCOFITS=self.getRawCOFITS()
        for eachSMFactor in self.smoothFactors:
            print "Smooting by factor of ",eachSMFactor
            self.smoothCubeSingle(rawCOFITS,  savePath=self.tmpPath,disFactor=eachSMFactor)


    def getFITSnamesByPrefix(self,path, prefix ):

        """
        Search all files according to
        :param path:
        :param prefix:
        :return:
        """
        fileList=[]

        for eachFactor in self.smoothFactors:


            searchStr = path+ prefix + "*{}*.fits".format( eachFactor )
            fileNames  = sorted(glob.glob(searchStr))  # sort this in names

            excludNoise=[]

            for eachFN in fileNames:
                if "noise" not in eachFN:
                    excludNoise.append(eachFN)

            if len(excludNoise)>1 :

                print "Multiple files found, wrong data ... exiting "
                sys.exit()
            fileList.append(   excludNoise[0] )


        return fileList


    def getTotalFlux(self,TBList):
        """
        return
        :param TBList:
        :return:
        """
        fluxList = []

        for eachTB in TBList:
            totalFlux=np.sum(eachTB["sum"])*self.getVelResolution() #convert to 0.2 km/s

            fluxList.append(totalFlux)

        return np.asarray(fluxList)

    def getTotalArea(self,TBList):
        """
        return
        :param TBList:
        :return:
        """
        areaList = []

        for eachTB in TBList:

            totalArea=np.sum(eachTB["area_exact"])  #convert to 0.2 km/s
            areaList.append(totalArea)

        return np.asarray(areaList)



    def getTBFixSmFactor(self,smFactor=1.0, dbParameters="S2P4Con1" ):
        """
        smFactor fixed, but change noiseFactor
        :param smFactor:
        :return:
        """

        tbList=[]
        noiseList= []
        for eachNF in self.noiseFactors:
            searchStr="{}*SmFactor*{}*noiseFactor*{}*{}*_Clean.fit".format(self.calCode,smFactor,eachNF,dbParameters)

            searchStr=os.path.join(self.tmpPath,searchStr)
            tbFile=glob.glob( searchStr )



            if len(tbFile) ==1:


                tbList.append(  Table.read( tbFile[0] ) )
                noiseList.append(   eachNF  )

            else:
                print "Multiple data found please check your data"
                print tbFile



        return np.asarray( tbList ) , np.asarray( noiseList )





    def produceAllRMSFITS(self):
        """
        produce RMS fits for all fits files
        :return:
        """
        for eachCode in self.allCodeList:
            fitsFile=self.getRawCOFITS(eachCode)
            saveName="{}_RMS.fits".format(eachCode)
            doFITS.getRMSFITS(fitsFile, saveName )

            if 0:#remove nan value
                rmsData,rmsHead=doFITS.readFITS(saveName)

                centerData=rmsData[:,200:2600]

                meanRMS=np.nanmean(centerData)

                rmsData[np.isnan(rmsData)]=meanRMS
                if "12" in eachCode: #less than 0.1 K is due to insufficent negative data
                    rmsData[rmsData<0.1]=meanRMS

                fits.writeto(saveName,rmsData,overwrite=True)

    def getMeanRMS(self):

        isRaw = self.isRaw( self.calCode )

        if "12" in self.calCode and not isRaw:
            return self.MWISPrmsCO12

        if "13" in self.calCode  and not isRaw:
            return self.MWISPrmsCO13

        if "18" in self.calCode  and not isRaw:
            return self.MWISPrmsCO18

        if "12" in self.calCode and  isRaw:
            return self.MWISPrmsRawCO12

        if "13" in self.calCode and  isRaw:
            return self.MWISPrmsRawCO13

        if "18" in self.calCode and  isRaw:
            return self.MWISPrmsRawCO18

        return None



    def cleanFITSsigma2(self,FITSfile,cutoff=2,minPts=4,contype=1  ):
        """

        #an updated version of DBSCAN code, which is more easily to use

        :param FITSfile:
        :param cutoff:
        :param minPts:
        :param contype:
        :return:
        """

        smFactor=self.getSmoothFactor(FITSfile)

        doMWdbscan = MWISPDBSCAN()
        doMWdbscan.rawCOFITS =  FITSfile
        doMWdbscan.rmsFITS = self.getRMSFITS()

        doMWdbscan.setDBSCANParameters( cutoff_sigma=cutoff,minPts=minPts,connectivity=contype)
        doMWdbscan.processPath = self.tmpPath

        doMWdbscan.computeDBSCAN()
        doMWdbscan.getCatFromLabelArray(doClean=True)

        #only produce clean fits for smFactor 1.0
        if smFactor==1.0:
            doMWdbscan.produceCleanFITS()



    def get3DRMSData(self,rmsData,Nz):

        """
        duplicate the rmsData
        :param rmsData:
        :param Nz:
        :return:
        """

        rmsDataMimic3d=np.asarray([rmsData])
        rmsData3D=np.repeat(rmsDataMimic3d,Nz,axis=0 )
        return rmsData3D

    def getMaskedRawCOFITS(self,maskSigma=2):

        rawCOFITS=self.getRawCOFITS()
        pathRaw, nameRaw =os.path.split(rawCOFITS)
        maskFITSName =  "downTo_{}sigma" + nameRaw

        maskedFITSFull=os.path.join(pathRaw,maskFITSName)

        if os.path.isfile(maskedFITSFull):
            return  maskedFITSFull

        else:
            return None


    def produceMaskSigma2(self,FITSFile,cutoffSigma=2):
        """
        this code, get raw fits, and rms fits, to mask the data to 2 sgima,level, which would be used to search molecular clouds
        fits would be removed after cleaning
        :param cutoffSigma:
        :return:
        """

        #rawCOFITS=self.getRawCOFITS()

        rmsFITS=self.getRMSFITS()

        print "Processing CO and FITS"
        print FITSFile
        print rmsFITS

        rawCOData,rawCOHead=doFITS.readFITS( FITSFile )

        rmsData,rmsHead = doFITS.readFITS(rmsFITS)
        Nz,Ny,Nx=rawCOData.shape

        rmsData3D= self.get3DRMSData(rmsData,Nz)

        
        #masked the ata
        rawCOData[rawCOData<cutoffSigma*rmsData3D]=np.nan

        #save the files

        saveName=self.tmpPath+"downTo_{}sigma".format(cutoffSigma)+os.path.basename(rawCOFITS)

        fits.writeto(saveName,rawCOData,header=rawCOHead,overwrite=True)

        return saveName

    def calRMSSingle(self,fitsName,onlyValue=False):
        """
        the input is a fits file
        :param filtsName:
        :return:
        """

        fitsPath,fitsNameSplit = os.path.split(fitsName)

        saveName= os.path.join(  fitsPath, "RMS_"+fitsNameSplit )


        if onlyValue:
            return doFITS.getRMSFITS(fitsName,"", returnRMSValue=True )

        else:
            doFITS.getRMSFITS(fitsName,saveName)



    def calRMSList(self,fitsFile ):
        """
        useThis instead of calRMS single
        :param fitsFileList:
        :return:
        """

        rmsList=[]

        #first check if the input fitsFile is
        if type(fitsFile) is list:

            for eachF in fitsFile:

                print "Processing ", eachF

                rmsOneFile = self.calRMSSingle(eachF,onlyValue=True)

                rmsList.append(rmsOneFile)

        return rmsList



    def smoothCubeSingle(self,FITSfile,  disFactor=2, savePath=None  ):
        """
        :param FITSfile:
        :param rawBeam:
        :param targetBeam:
        :return:
        """

        rawBeamSize = self.getBeamSize() # 52. / 60,
        fitsNameBase= os.path.basename(FITSfile)

        namePre,nameSuf=os.path.splitext(fitsNameBase)

        suffix=  '_SmFactor_{}'.format(disFactor)

        saveName=self.calCode+namePre+suffix+nameSuf # #self.addFixToFITSName(FITSfile, suffix,  prefix=False)
        #saveName= #self.addFixToFITSName(saveName, self.calCode,  prefix=True)


        ## ##  #

        if savePath!=None:
            saveName = savePath+ saveName
        else:
            saveName = self.tmpPath+ saveName


        if disFactor==1.0: #if the factor is 1, there is no need to calculate
            data,head=doFITS.readFITS(FITSfile)
            fits.writeto(saveName,data,header=head,overwrite = True )

            return


        targetBeam=disFactor*rawBeamSize

        processCube = SpectralCube.read( FITSfile )
        processCube.allow_huge_operations = True

        RawBeam = radio_beam.Beam(major=rawBeamSize * u.arcmin, minor=rawBeamSize * u.arcmin, pa=0 * u.deg)
        processCube.beam =  RawBeam

        beamTarget = radio_beam.Beam(major=targetBeam * u.arcmin, minor= targetBeam * u.arcmin, pa=0 * u.deg)


        new_cube = processCube.convolve_to(beamTarget)

        new_cube.write(saveName,overwrite=True)

        doFITS.converto32bit(saveName,saveName)


    def addFixToFITSName(self, fitsName , addNote, prefix=False ):
        """
        Insert prefix to to name of the fitsName
        :param fitsName:
        :return:
        """


        #first split the name


        path,rawName=os.path.split( fitsName )

        rawNamePre,rawNameExt=os.path.splitext(rawName)


        if prefix:
            hybirdName =  addNote+rawNamePre+rawNameExt

        else:
            hybirdName = rawNamePre+addNote+rawNameExt


        return os.path.join(path,hybirdName)



    def addNoseToAllSmFiles(self,doClean=False):
        """
        search all smoothed files,
        :return:
        """

        #step 1, get all smoothed files

        fitsList = self.getFITSnamesByPrefix(self.tmpPath,self.calCode)

        for eachFile in fitsList:
            for eachNFactor in self.noiseFactors:


                print eachNFactor,eachFile

                noiseFITS = self.addNoiseSingle(eachFile,eachNFactor)
                print noiseFITS,"????????????????????"

                if doClean:
                    self.cleanFITSsigma2(noiseFITS,removeFITS=True)
                    os.remove(noiseFITS) #if do clean, then remove the noisefits, only keep the table





    def addNoiseSingle(self,smoothedFITS,  noiseFactor=2, absoluteNoise=None , returnData=False, redo=False ):
        """

        :param smoothedFITS:
        :param noiseFactor:
        :return:
        """

        #this code add noise by a mean std obsoleted,
        #
        aaaaa

        saveSuffix= "_noiseFactor_{}".format( float(noiseFactor) )

        saveName=self.addFixToFITSName(smoothedFITS,saveSuffix,prefix=False)

        if not redo and os.path.isfile( saveName ):
            return




        if noiseFactor==1.0:
            data, head = myFITS.readFITS(smoothedFITS)

            fits.writeto(saveName, data, header=head, overwrite=True)

            return saveName


        smoothedRMS=doFITS.getRMSFITS(smoothedFITS, ""  , returnRMSValue=True )

        if absoluteNoise==None:
            targedRMS=noiseFactor*smoothedRMS
        else:
            targedRMS= absoluteNoise
            saveSuffix = "_noise_{}absK".format(float(absoluteNoise))

            saveName=self.addFixToFITSName(smoothedFITS,saveSuffix,prefix=False)


        print "the smoothed rms is {}, the targed rms is {}".format( smoothedRMS,  targedRMS)


        convolvedStd=np.sqrt(targedRMS**2- smoothedRMS**2 )
        print "The convolved sts is {} K ".format( convolvedStd )

        data,head=myFITS.readFITS(  smoothedFITS )

        noise=np.random.normal(0, convolvedStd,data.shape)

        observedData=data+noise

        #del head["CELLSCAL"]
        #del head["BSCALE"]
        #del head["BZERO"]
        if returnData:

            return observedData



        fits.writeto( saveName , observedData,header= head,overwrite=True )
        doFITS.converto32bit(saveName,saveName)

        return saveName

    def getSmFITSFileList(self,calCode=None):

        fileList=[]

        processCode=self.calCode

        if calCode!=None:
            processCode = calCode

        for i in self.smoothFactors:

            searchStr =  "{}*SmFactor_{:.1f}.fits".format(processCode, i  )
            searchStr = os.path.join( self.tmpPath, searchStr)

            #print searchStr

            fitsName=glob.glob(searchStr)

            if len(fitsName )==1:
                fileList.append(fitsName[0] )

            else:
                print fitsName
                print "{} fits found (none or multiple, reject this search)...".format(len(fitsName))


        return  fileList

    def getSmFITSFileSingle(self, smFactor , calCode=None):

        processCode=self.calCode

        if calCode!=None:
            processCode = calCode


        searchStr =  "{}*{}{:.1f}.fits".format(processCode,self.smoothTag, smFactor  )
        searchStr = os.path.join( self.tmpPath, searchStr)

        #print searchStr

        fitsName=glob.glob(searchStr)

        if len(fitsName )==1:
            #fileList.append(fitsName[0] )
            return fitsName[0]
        else:
            print fitsName
            print "{} fits found (none or multiple, reject this search)...".format(len(fitsName))

            return None

    def getProcessCod(self,calCode):
        processCode=self.calCode
        if calCode!=None:
            processCode = calCode

        return processCode

    def getSmoothListFixNoise(self, noiseFactor=0.0, calCode=None,getCleanTBFile=False,getCleanFITS=False):
        """
        fix noise factors, get a list of fits file accordin the smoooth factors

        :param noiseFactor:
        :param calCode:
        :return:
        """

        fileList=[]

        for eachSmFactor in self.smoothFactors:

            file=self.getSmoothAndNoiseFITSSingle(smFactor=eachSmFactor,noiseFactor=noiseFactor,calCode=calCode , getCleanTBFile=getCleanTBFile, getCleanFITS=getCleanFITS )

            fileList.append(file)

        return fileList



    def getSmoothListFixSmooth(self, smoothFactor=1.0, calCode=None,getCleanTBFile=False,getCleanFITS=False):
        """
        fix noise factors, get a list of fits file accordin the smoooth factors

        :param noiseFactor:
        :param calCode:
        :return:
        """

        fileList=[]

        for eachNoiseFactor in self.noiseFactors:

            file=self.getSmoothAndNoiseFITSSingle(smFactor=smoothFactor,noiseFactor=eachNoiseFactor,calCode=calCode,getCleanTBFile=getCleanTBFile, getCleanFITS=getCleanFITS )

            fileList.append(file)

        return fileList



    def getSmoothAndNoiseFITSSingle( self, smFactor=1.0, noiseFactor=0.0, dbscanCode="dbscanS2P4Con1" ,calCode=None, getCleanTBFile=False,getCleanFITS=False, getRawDBSCANFITS=False):
        """

        :param smFactor:
        :param noise:
        :return:
        """

        processCode=self.getProcessCod(calCode)
        searchStr = "{}*{}{:.1f}*{}{:.1f}.fits".format(processCode, self.smoothTag, smFactor, self.noiseTag,
                                                       noiseFactor)

        if getCleanTBFile:
            searchStr =  "{}*{}{:.1f}*{}{:.1f}{}_Clean.fit".format(processCode,self.smoothTag, smFactor ,self.noiseTag,noiseFactor,dbscanCode )

        if getCleanFITS:
            searchStr =  "{}*{}{:.1f}*{}{:.1f}{}_Clean.fits".format(processCode,self.smoothTag, smFactor ,self.noiseTag,noiseFactor,dbscanCode )

        if getRawDBSCANFITS:
            searchStr = "{}*{}{:.1f}*{}{:.1f}{}.fits".format(processCode, self.smoothTag, smFactor, self.noiseTag,  noiseFactor, dbscanCode)



        searchStr = os.path.join( self.tmpPath, searchStr)


        fitsName=glob.glob(searchStr)

        if len(fitsName )==1:
            #fileList.append(fitsName[0] )
            return fitsName[0]
        else:
            print fitsName
            print "{} fits found (none or multiple, reject this search)...".format(len(fitsName))

            return None



    def drawBeamEffect(self):
        """
        draw the effect of the sensitivity, use color to code different smooth factors
        :return:
        """


        fig = plt.figure(figsize=(10, 8))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]


        axBeam = fig.add_subplot(1,1, 1)


        #draw with color
        for eachNoiseFactor in self.noiseFactors:
            self.drawBeamEffectSingle(axBeam, noiseFactor= eachNoiseFactor )
        #self.drawBeamEffectSingle(axBeam, noiseFactor=1.5  )


        #self.drawNoiseEffectSingle(axSens,1.5)

        axBeam.set_ylabel(r"Flux ($\rm K\ km\ s$$^{-1}$ $\mathrm{\Omega}_\mathrm{A}$)")
        axBeam.set_xlabel(r"Beam Size (arcmin)")

        self.drawColorBarNoiseFactor(axBeam)

        fig.tight_layout()
        plt.savefig("effectBeam.pdf", bbox_inches='tight')
        plt.savefig("effectBeam.png", bbox_inches='tight', dpi=600)


    def drawBeamEffectSingle(self,ax,noiseFactor):

        tbList,beamFactorArray=self.getTBFixNoiseFactor(  noiseFactor=  noiseFactor )

        beamSizeArray= beamFactorArray*self.rawBeamSize


        #for each



        totalFlux=self.getTotalFlux(tbList)


        #better convert noise factor to Kelvin
        #get the noise factor

        #fitsRawFile=self.getFITSBySF(SF=smFactor   )

        #factor1RMS=self.calRMSSingle(fitsRawFile, onlyValue=True)
        #print "The rms for factor 1 is, ", factor1RMS

        plotColor = self.noiseColor.to_rgba(noiseFactor)

        ax.plot( beamSizeArray , totalFlux  ,'o-',   color=  plotColor  ,markersize=2 ,lw=1)




    def drawNoiseEffect(self):
        """
        draw the effect of the sensitivity, use color to code different smooth factors
        :return:
        """


        fig = plt.figure(figsize=(10, 8))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]


        axSens = fig.add_subplot(1,1, 1)


        #draw with color
        for eachFactor in  self.smoothFactors:
            self.drawNoiseEffectSingle(axSens, eachFactor )



        axSens.set_ylabel(r"Flux ($\rm K\ km\ s$$^{-1}$ $\mathrm{\Omega}_\mathrm{A}$)")
        axSens.set_xlabel(r"noise RMS (K)")

        self.drawColorBarSmooth(axSens)

        fig.tight_layout()
        plt.savefig("effectSensitivity.pdf", bbox_inches='tight')
        plt.savefig("effectSensitivity.png", bbox_inches='tight', dpi=600)



    def drawColorBarSmooth(self,axCon3):

        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(axCon3)
        cax1 = divider.append_axes("right", size="3%", pad=0.05)
        self.smoothColor._A=[]
        cb = plt.colorbar( self.smoothColor , cax=cax1)

        cax1.set_ylabel("Beam Size (arcmin)")

    def drawColorBarNoiseFactor(self,axCon3):

        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(axCon3)
        cax1 = divider.append_axes("right", size="3%", pad=0.05)
        self.noiseColor._A=[]
        cb = plt.colorbar( self.noiseColor , cax=cax1)

        cax1.set_ylabel("Noise RMS increase factor")




    def removeNoiseFiles(self,fitsNameList):
        """

        :param fitsNameList:
        :return:
        """

        returnList=[]

        for eachF in fitsNameList:
            if self.noiseStr in eachF:
                continue
            else:
                returnList.append(eachF)

        return returnList


    def getFileByStr(self,path,nameStr):


        searchStr= os.path.join( path, nameStr )
        foundFiles= glob.glob(searchStr)

        if len( foundFiles )==1:
            return foundFiles[0]

        else:
            print searchStr
            print "Please check you data!"
            return None


    def getAbsNoiseFITSName(self,absK=0.5):

        fitsList= []

        for eachSF in self.smoothFactors:
            fileName="{}*SmFactor_{}*{}absK.fits".format(self.calCode, eachSF ,absK )

            fitsList.append( self.getFileByStr(self.tmpPath,fileName) )

        return fitsList



    def getAbsNoiseCleanTBList(self,absK=0.5, dbscanCode="dbscanS2P4Con1"):

        tbList= []
        tbFileList= []
        for eachSF in self.smoothFactors:
            fileName="{}*SmFactor_{}*{:.2f}absK{}_Clean.fit".format(self.calCode, eachSF ,absK, dbscanCode )

            TBName= self.getFileByStr(self.tmpPath,fileName)

            if TBName!=None:

                tbList.append(Table.read(TBName) )

            else:

                tbList.append(TBName )
            tbFileList.append(TBName)
        return tbList,tbFileList


    def sumVoxels(self,fitsName,downToK=1.0):
        """
        The unite of downToK is
        :param downToK:
        :return:
        """

        if type(fitsName)==list:
            #
            sumList=[]

            for eachF in  fitsName:
                print "Processing ",eachF
                data, head = myFITS.readFITS(eachF)
                self.maskNoiseEdge(data)

                sumV= np.nansum(data[data >= downToK])*self.getVelResolution()

                sumList.append( sumV )

            return np.asarray( sumList )


        else:
            data,head=myFITS.readFITS(fitsName)
            self.maskNoiseEdge(data)
            #print np.nansum(data)
            return np.nansum(data[data>=downToK])*self.getVelResolution()


    def maskNoiseEdge(self,data):
        """
        Only apply to local molecular clouds
        :param data:
        :return:
        """

        data[ :, 1063:  ,0:56 ,]=0
        data[:, 1003: , 2815: ]=0


    def compareSensitivity(self, smoothedFITS ):

        """
        This function is used test different methods of cutoff
        the first method
        :param smoothedFITS:
        :return:
        """
        #compare direct cutoff and the cutoff after increasing the

        data,head= doFITS.readFITS( smoothedFITS )

        #cutOffK=[2,3,4,5,6,7,8,9,10] #in K
        cutOffK=[ 10 ] #in K

        for eachK in cutOffK:

            dirrectCutValue = np.nansum(data[data>eachK])*self.getVelResolution()

            #if we take the dirrect cutValue as the 2 sigma
            noiseData=self.addNoiseSingle( smoothedFITS ,   absoluteNoise= eachK/4. , returnData=True)

            noiseDataSelect = noiseData[noiseData>eachK]

            print "Total Number of summed voxels, ", len( noiseDataSelect )

            noiseCutValue = np.nansum(noiseDataSelect)*self.getVelResolution()

            diff=   noiseCutValue - dirrectCutValue
            print  dirrectCutValue,  noiseCutValue,   diff/dirrectCutValue


    def fluxAboveCutoff(self,data,cutoff):
        """
        sum all flux above the cutoff threshold
        :param data:
        :param cutoff:
        :return:
        """

        return np.nansum(data[data>cutoff])*self.getVelResolution()



    def testEqualCutOff(self):
        """
        Find the equivalent cut off, that would yield a equal flux with the cloud flux
        :return:
        """

        #for certain test code , find equivalent cut off sigma value, and use  the average value as the cut off,
        # is this cut off vlaue diffierent, for different region and lines?

        ##find all sm files

        ##find all TB files

        #smFiles=self.getSmFITSFileList()
        smFiles = self.getSmoothListFixNoise( )

        ceanTBFiless = self.getSmoothListFixNoise(getCleanTBFile=True)

        eqSigmaList = []

        for i in range(len(smFiles)):
            fitsName = smFiles[i]
            tbFileName = ceanTBFiless[i]

            print fitsName
            print tbFileName

            eqDownToSigma = self.getEqualCutOff(fitsName, tbFileName) / self.getMeanRMS()
            print eqDownToSigma
            eqSigmaList.append(eqDownToSigma)

        print eqSigmaList
        print  "the equal sigma list is ", np.mean(eqSigmaList)

    def getEqualCutOff(self,fitsFile,TBFile):

        """
        get the equivalent cutoff of fitsFile to get the equivalent flux with DBSCAN catalg
        :param fitsFile:
        :param TBFile:
        :return:
        """


        tbFile=Table.read(TBFile)

        dbscanFlux = np.sum(tbFile["sum"])*self.getVelResolution()

        precision=0.0001 #percent



        downToKUp = 3. #upper, the lower is 3 K
        downToKLow = 1.0 #upper, the lower is 1 K


        data, head = myFITS.readFITS(fitsFile)
        #self.maskNoiseEdge(data)



        while 1:
            #calcualate cutoff rms
            downToKTry = (downToKUp + downToKLow) / 2.
            #print "Trying downToK, ", downToKTry
            sumV = np.nansum(data[data >= downToKTry]) *self.getVelResolution()

            ###
            diff= abs( sumV-dbscanFlux)/dbscanFlux
            #print diff
            if diff<=precision:
                break

            else:
                if sumV > dbscanFlux: # the sumV is too large , incrase the cutoff
                    downToKLow= downToKTry

                else: # sumV is too small, decrease the try value by dowong
                    downToKUp= downToKTry
            print downToKTry

        return downToKTry

    def convertFunction(self,x):

        return np.sqrt(x)

    def guessNoiseByCode(self):

        noise = self.MWISPrmsCO12

        if "CO13" in self.calCode:
            noise=self.MWISPrmsCO13

        if "CO18" in self.calCode:
            noise=self.MWISPrmsCO18

        return noise

    def getNoiseDataList( self ):

        targetNoise=self.guessNoiseByCode()

        noiseFiles=self.getAbsNoiseFITSName( targetNoise )

        dataList=[]


        for eachF in noiseFiles:

            dataN,headN=doFITS.readFITS( eachF )

            dataList.append(dataN)

        return dataList


    def fittingFFAndSize(self,useFunction,sizeArray,ffArray,ffError, sizeError=None,useODR=False):
        """
        test the number of parameters of userFUnction, to determine which function should be used
        :param x:
        :param y:
        :return:
        """


        if not useODR:


            parametersN=useFunction.func_code.co_argcount
            print parametersN,"?????????????"
            if parametersN== 4:
                params, paramas_covariance = optimize.curve_fit(useFunction, sizeArray, ffArray, sigma=ffError,   absolute_sigma=True, p0=[1, 0.5, 1])

            if parametersN == 3:
                #params, paramas_covariance = optimize.curve_fit(useFunction, sizeArray, ffArray, sigma=ffError,    absolute_sigma=True, p0=[  0.5, 1])
                params, paramas_covariance = optimize.curve_fit(useFunction, sizeArray, ffArray, sigma=ffError,    absolute_sigma=True, p0=[  1, 0.5])


            if parametersN == 2:
                params, paramas_covariance = optimize.curve_fit(useFunction, sizeArray, ffArray, sigma=ffError,    absolute_sigma=True, p0=[ 4])



        #fitting a function that describe the
            errors = np.sqrt(np.diag(paramas_covariance))
            return params,errors


        if useODR:
            #fitting function with ODR

            ffModel = odrpack.Model(testffAndSizeFuncODR)
            if sizeError is None:
                mydata = odrpack.RealData(sizeArray, ffArray, sy=ffError)
            else:
                mydata = odrpack.RealData(sizeArray , ffArray, sx=sizeError, sy=ffError)


            myodr = odrpack.ODR(mydata, ffModel, beta0=[4])
            myOutPut = myodr.run()
            # myOutPut.pprint()

            return myOutPut.beta, myOutPut.sd_beta

        return 0,   0



    def fittingFFSizeRelationOdr(self,sizeArray,ffArray,sizeError, ffError):
        """
        #for over all fiting

        :param sizeArray:
        :param ffArray:
        :param ffError:
        :return:
        """

        #pass
        print "Fitting filling factor and size relationship"
        ffModel = odrpack.Model(ffAndSizeFunc  )

        mydata = odrpack.RealData(sizeArray,ffArray,  sx= sizeError, sy=ffError)

        myodr = odrpack.ODR(mydata, ffModel, beta0=[ 0.5, 1.])
        myOutPut = myodr.run()
        myOutPut.pprint()

        return   myOutPut.beta,  myOutPut.sd_beta



    def ffOdr(self,x,y,yError,xError=None):
        """
        fitting
        :param x:
        :param y:
        :param yError:
        :return:
        """

        ffModel = odrpack.Model(ffFunctionOdr  )
        if xError is None:
            mydata = odrpack.RealData(x,y,   sy=yError)
        else:
            mydata = odrpack.RealData(x,y,  sx=xError, sy=yError)

        averageY=np.mean(y)
        myodr = odrpack.ODR(mydata, ffModel, beta0=[  averageY , 0.5, averageY ] )
        myOutPut = myodr.run()
        #myOutPut.pprint()

        return   myOutPut.beta,  myOutPut.sd_beta




    def getFillingFactorAndDraw(self,beamList,fluxList,fluxError,calID,drawFigure=False,saveTag=""):
        """

        :param beamList:
        :param fluxList:
        :return:
        """

        x=np.asarray( beamList )

        y = np.asarray( fluxList )

        #print x
        #print y
        #try:



        if fluxError is not None:


            params, paramas_covariance = optimize.curve_fit(ffFunction, x, y, sigma=fluxError, absolute_sigma=True,  p0=[np.mean(y), 0.5, np.mean(y)])

        else:
            params, paramas_covariance = optimize.curve_fit(ffFunction, x, y,   p0=[np.mean(y), 0.5, np.mean(y)])


        errors = np.sqrt(np.diag(paramas_covariance))

        #params,errors=self.ffOdr(x,y,fluxError) #dot no use Odr for filing factor fitting, fore there is no error in axis


        #except:

            #return 0, 0, [[0,0,0],[0,0,0]]



        #errors = np.sqrt(np.diag(paramas_covariance))

        fittingParaAndError = [params, errors]

        # for i in range(len(params)):
        # print params[i], errors[i],"Error in percentage", errors[i]/params[i]
        cfaBeam=8.5
        cfaFilling =        ffFunction(cfaBeam, params[0], params[1], params[2]) / ffFunction(0, params[0], params[1], params[2])
        wmsipFilling = ffFunction(self.getBeamSize( ), params[0], params[1], params[2]) / ffFunction(0, params[0], params[1],   params[2])
        if not drawFigure:
            return wmsipFilling, cfaFilling, fittingParaAndError

        #drawTheFigure
        fig = plt.figure(figsize=(10, 8))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})

        if 0:
            mpl.rcParams['text.latex.preamble'] = [
                r'\usepackage{tgheros}',  # helvetica font
                r'\usepackage{sansmath}',  # math-font matching  helvetica
                r'\sansmath'  # actually tell tex to use it!
                r'\usepackage{siunitx}',  # micro symbols
                r'\sisetup{detect-all}',  # force siunitx to use the fonts
            ]

        axFitting = fig.add_subplot(1,1, 1)
        #axFitting.scatter(x, y, s=15, color='red', label="The noise RMS is 0.5 K")
        #axFitting.scatter(x, x*y, s=15, color='red', label="The noise RMS is 0.5 K")
        #axFitting.scatter(  x  ,  y , s=15,  color='red'   ) #
        axFitting.errorbar(x,   y , yerr= fluxError ,    markersize=2 , fmt='o', c= 'red' ,     capsize=0.5,  elinewidth=1 , lw=1 )

        fittingX = np.arange(0, np.max(x), 0.01)
        # axFitting.plot( fittingX  ,  ffFunction(fittingX,params[0], params[1] , params[2]   ), color='blue'  )
        axFitting.plot(fittingX, ffFunction(fittingX, *params), color='blue', lw=1.0)

        fillingWMISP = "The filling factor ({}) of MWISP is: {:.3f}".format(self.calCode, wmsipFilling)
        fillingCfA = "The filling factor ({}) of CfA is: {:.3f}".format(self.calCode, cfaFilling)

        #three paramters

        paramterValues="a: {:.2f}(+/-{:.2f}), b: {:.2f}(+/-{:.2f}), c: {:.2f}(+/-{:.2f})".format(params[0],errors[0],params[1],errors[1],params[2],errors[2])

        # at = AnchoredText(r"The noise rms is {} K (1$\sigma$)".format(self.MWISPrmsCO12), loc=1, frameon=False)
        # axFitting.add_artist(at)

        # if cfaBeam> np.max(x):
        # cfaFilling= 0

        at = AnchoredText(fillingWMISP + "\n" + fillingCfA+"\n"+paramterValues, loc=1, frameon=False)
        axFitting.add_artist(at)
        ###

        axFitting.axvline(x=0, ls="--", color='black')

        # draw CfA resolution
        axFitting.axvline(x=cfaBeam, ymax=0.2, ls="--", color='green', label="CfA resolutuion (8.5 arcmin)")
        axFitting.legend(loc=5)
        ###########################
        saveTag= "{}_factorFitting_ID{}{}".format(self.calCode,calID,saveTag)
        if calID!=0:
            plt.savefig(self.figurePath+"{}.png".format( saveTag ), bbox_inches='tight', dpi=200)
            
        else:
            plt.savefig(self.paperFigurePath+"{}.png".format( saveTag ), bbox_inches='tight', dpi=200)
            plt.savefig(self.paperFigurePath+"{}.pdf".format( saveTag ), bbox_inches='tight' )

            

        plt.close(fig )
        gc.collect()



        return wmsipFilling, cfaFilling,  fittingParaAndError

    def getFluxAndErrorList(self,TBList):
        """
        TBList, shold be an tb file name list
        :param TBList:
        :return:
        """

        fluxList=[]
        errorList = []


        velReslution=self.getVelResolution()
        meanRMS=self.getMeanRMS() #*velReslution

        for eachTB in TBList:


            processTB=Table.read(eachTB)
            #remove wrong clouds caused by bad channels

            processTB=self.removeFakeClouds(processTB)



            fluxList.append( np.sum(processTB["sum"])*velReslution )
            print eachTB, np.sum(processTB["sum"])*velReslution
            fluxError = np.sqrt(np.sum(processTB["pixN"]))*meanRMS*velReslution

            errorList.append(fluxError)

        return np.asarray(fluxList), np.asarray(errorList)




    def callFillingFactorAllSM(self   ):

        """
        Get clean TB files to draw total flux change, getFillingFactors
        :param calCode:
        :param absK:
        :return:
        """


        mwispFillingList=[]
        cfaFillingList=[]
        fittingParaList=[]



        for eachCode  in  self.allRawCodeList:

            self.calCode=eachCode

            beamSizeArray = self.getBeamSize() * self.smoothFactors

            smTBFiles = self.getSmoothListFixNoise(getCleanTBFile=True)

            fluxList, fluxError = self.getFluxAndErrorList(smTBFiles)




            mwispFilling, cfaFilling, fittingParaAndError = self.getFillingFactorAndDraw(beamSizeArray, fluxList,fluxError, 0, drawFigure=True)

            mwispFillingList.append( mwispFilling )
            cfaFillingList.append( cfaFilling )
            fittingParaList.append( fittingParaAndError )

            x = self.getBeamSize() * self.smoothFactors
            np.save(self.calCode + "X", x)
            np.save(self.calCode + "Y", fluxList)



        np.save("mwispFillingRaw", mwispFillingList )
        np.save("cfaFillingRaw", cfaFillingList )
        np.save("fittingParaRaw", fittingParaList )


        return

        aaaaaaaaaaaaa


        ##

        TBCleanTB,TBFiles=self.getAbsNoiseCleanTBLIst(absK=absK)
        print TBFiles[0]
        y= self.getTotalFlux( TBCleanTB )

        saveTag= "{}_factorFitting_flux".format(self.calCode)

        if useArea:
            y= self.getTotalArea( TBCleanTB )
            saveTag= "{}_factorFitting_area".format(self.calCode)


        x =self.rawBeamSize*self.smoothFactors
        np.save(self.calCode+"X",x)
        np.save(self.calCode+"Y",y)

        #x=np.exp(x)

        ################################################
        #select=y>0
        #x=x[select]
        #y=y[select]

        #x= self.convertFunction( x )
        #x= x**2

        #y= np.log(y)  #self.getTotalFlux( TBCleanTB )

        ########
        fig = plt.figure(figsize=(10, 8))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})

        if 0:

            mpl.rcParams['text.latex.preamble'] = [
                r'\usepackage{tgheros}',  # helvetica font
                r'\usepackage{sansmath}',  # math-font matching  helvetica
                r'\sansmath'  # actually tell tex to use it!
                r'\usepackage{siunitx}',  # micro symbols
                r'\sisetup{detect-all}',  # force siunitx to use the fonts
            ]

        axFitting = fig.add_subplot(1,1, 1)
        #axFitting.scatter(x, y, s=15, color='red', label="The noise RMS is 0.5 K")
        #axFitting.scatter(x, x*y, s=15, color='red', label="The noise RMS is 0.5 K")
        axFitting.scatter(  x  ,  y , s=15,  color='red', label="The noise RMS is {:.2f} K".format( absK ) )


        if 0: #curving fitting


            params, paramas_covariance=optimize.curve_fit(ffFunction, x,y, p0=[10000, 0.5,   np.mean(y)  ] )
            print params
            fittingX=np.arange(0, np.max(x),0.01)
            axFitting.plot( fittingX  ,  ffFunction(fittingX,params[0], params[1] , params[2]   ), color='blue'  )

        if 1: #do not draw xy

            cfaFilling = 1
            wmsipFilling = 1
            cfaBeam= 8.5 #arcmin
            if 0:
                pnDegree = 5

                z1 = np.polyfit(x, y, pnDegree )
                p1 = np.poly1d(z1)

                cfaFilling= p1(  cfaBeam )/p1(0)


                fittingX=np.arange(-0.2, np.max(x),0.1)
                axFitting.plot( fittingX,  p1(fittingX) ,color='blue' ,lw=1 ,label="Polynomical fitting with degree {}".format(pnDegree ) )

                fillingWMISP =  "The filling factor ({}) of MWISP is: {:.3f}".format(self.calCode, p1(self.rawBeamSize)/p1(0) )
                fillingCfA = "The filling factor ({}) of CfA is: {:.3f}".format(self.calCode, cfaFilling )


            else:
                params, paramas_covariance=optimize.curve_fit(ffFunction, x,y, p0=[10000, 0.5,   np.mean(y)  ] )
                errors=np.sqrt(np.diag(paramas_covariance))

                fittingParaAndError=[ params,  errors  ]

                #for i in range(len(params)):
                    #print params[i], errors[i],"Error in percentage", errors[i]/params[i]

                cfaFilling =    ffFunction( cfaBeam , params[0], params[1] , params[2]   ) /  ffFunction(0, params[0], params[1] , params[2]   )
                wmsipFilling =   ffFunction(self.rawBeamSize,  params[0], params[1] , params[2]    ) /  ffFunction(0, params[0], params[1] , params[2]   )

                fittingX=np.arange(0, np.max(x),0.01)
                #axFitting.plot( fittingX  ,  ffFunction(fittingX,params[0], params[1] , params[2]   ), color='blue'  )
                axFitting.plot( fittingX  ,  ffFunction(fittingX,*params   ), color='blue' ,lw=1.5 )

                fillingWMISP = "The filling factor ({}) of MWISP is: {:.3f}".format(self.calCode,   wmsipFilling )
                fillingCfA = "The filling factor ({}) of CfA is: {:.3f}".format(self.calCode, cfaFilling)

                #at = AnchoredText(r"The noise rms is {} K (1$\sigma$)".format(self.MWISPrmsCO12), loc=1, frameon=False)
                #axFitting.add_artist(at)


            #if cfaBeam> np.max(x):
                #cfaFilling= 0


            at = AnchoredText( fillingWMISP+"\n" +  fillingCfA , loc=1, frameon=False)
            axFitting.add_artist(at)

            axFitting.axvline(  x=0  ,  ls="--" ,  color='black'   )


            #draw CfA resolution
            axFitting.axvline(  x=cfaBeam  , ymax= 0.2, ls="--" ,  color='green' ,  label="CfA resolutuion (8.5 arcmin)" )
            axFitting.legend(loc= 5 )

        ######
        if useArea:
            axFitting.set_ylabel(r"Angular area (square arcmin)")

        else:
            axFitting.set_ylabel(r"Flux ($\rm K\ km\ s$$^{-1}$ $\mathrm{\Omega}_\mathrm{A}$)")


        axFitting.set_xlabel(r"Beam Size (arcmin)")



        fig.tight_layout()


        #axFitting.set_yscale('log')
        #axFitting.set_xscale('log')

        plt.savefig(self.ffFigurePath+"{}.pdf".format( saveTag  ), bbox_inches='tight')
        plt.savefig(self.ffffFigurePath+"{}.png".format( saveTag ), bbox_inches='tight', dpi=600)
        #the flux is saved in npy
        return   wmsipFilling, cfaFilling  , fittingParaAndError


    def printFillingCat(self):
        """
        print the filling factor for each arm, each lines and together with the fitting parameters
        :return:
        """

        ###

        mwispFillingList=np.load("mwispFillingRaw.npy")
        cfaFillingList=np.load("cfaFillingRaw.npy")
        paraFitting = np.load("fittingParaRaw.npy")


        for i in range(len(mwispFillingList)) :



            armStr= self.allArmList[i]

            lineStr= self.allLineList[i]


            factorStr= "{:.3f}".format( mwispFillingList[i] )
            para,paraError = paraFitting[i]

            aStr= self.getScitificStr(para[0], paraError[0])
            bStr=  "{:.3f} $\pm$ {:.3f}".format( para[1],paraError[1]  )
            #cStr= "{:.2f} $\pm$ {:.2f}".format( para[2],paraError[2]  )


            if para[2]<10:
                cStr= "{:.3f} $\pm$ {:.3f}".format( para[2],paraError[2]  )
            else:
                cStr=  self.getScitificStr(para[2], paraError[2])


            #print lineStr,armStr, "?????????????????????????"

            if lineStr==self.lineStrList[-1] and armStr==self.allArmList[-1]:


                lintStr="   & {} & -- & -- & -- & --  \\\\".format(   lineStr    )


            else:

                if i % 3 != 1:
                    armStr = ""

                lintStr="{} & {} & {} & {} & {} & {}  \\\\".format( armStr , lineStr , factorStr,  aStr, bStr, cStr   )



            print lintStr

            if i % 3 == 2 and i< len( mwispFillingList)-1:
                print "\\hline"

    def getScitificStr(self,value,valueError):
        """

        :param value:
        :param valueError:
        :return:
        """

        sciNote= "{:.3e}".format(value)
        sciPart=sciNote[5:]

        order = float( "1.0"+sciPart)
        orderN= int( np.log10( order) )

        valueFloat = value / order
        errorFloat= valueError/order


        finalStr= "({:.3f}$\\pm${:.3f})$\\times10^{{{}}}$".format(valueFloat ,errorFloat, orderN)
        return finalStr





    def printFluxCat(self):
        """

        :return:
        """
        mwispFillingList=np.load("mwispFillingRaw.npy")
        cfaFillingList=np.load("cfaFillingRaw.npy")
        paraFitting = np.load("fittingParaRaw.npy")





        for i in range(len(mwispFillingList)) :



            armStr= self.allArmList[i]

            lineStr= self.allLineList[i]


            factorStr= "{:.3f}".format( mwispFillingList[i] )
            para,paraError = paraFitting[i]

            aStr= self.getScitificStr(para[0], paraError[0])
            bStr=  "{:.3f} $\pm$ {:.3f}".format( para[1],paraError[1]  )
            #cStr= "{:.2f} $\pm$ {:.2f}".format( para[2],paraError[2]  )


            if para[2]<10:
                cStr= "{:.3f} $\pm$ {:.3f}".format( para[2],paraError[2]  )

            else:
                cStr=  self.getScitificStr(para[2], paraError[2])

            if i % 3 != 1:
                armStr=""

            codeI=self.allRawCodeList[i]
            X=np.load(codeI+"X.npy")
            Y=np.load(codeI+"Y.npy")

            strFlux = self.getLatexFromArray(Y)

            lintStr="{} & {} & {}  \\\\".format( armStr , lineStr , strFlux  )



            print lintStr

            if i % 3 == 2 and i< len( mwispFillingList)-1:
                print "\\hline"



    def getLatexFromArray(self , valueList ):
        """

        :param valueList:
        :return:
        """
        str=""
        for i in range( len(valueList)  ):
            value=valueList[i]/1000.
            latextFlux= self.latexFlux(value)
            if i==0:

                str=str +"{}".format( latextFlux )

            else:

                str=str +" & {}".format( latextFlux )

        return str


    def latexFlux(self,value):


        if value>100:
            return "{}".format(int(value))


        if value>10:
            return "{:.1f}".format( value )

        if value>1:
            return "{:.2f}".format( value )

        return "{:.3f}".format( value )



    def getScitificLatexSingle(self,value):
        """

        :param value:
        :return:
        """
        sciNote= "{:.3e}".format(value)
        sciPart=sciNote[5:]

        order = float( "1.0"+sciPart)
        orderN= int( np.log10( order) )


        valueFloat = value / order

        finalStr= "{:.2f}$\\times10^{{{}}}$".format(valueFloat ,  orderN)

        return finalStr




    def recordFillingFactorData(self):
        """
        print the filling factor of
        :return:
        """

        mwispFillingList=[]
        cfaFillingList=[]
        fittingParaList=[]



        for eachCode, eachNoise in zip(self.allCodeList,self.allNoiseList):

            self.calCode=eachCode
            mwispFilling ,cfaFilling, fittingParaAndError= self.drawFillingFactor(absK=eachNoise, useArea=False )

            mwispFillingList.append( mwispFilling )
            cfaFillingList.append( cfaFilling )
            fittingParaList.append( fittingParaAndError )

        np.save("mwispFilling", mwispFillingList )
        np.save("cfaFilling", cfaFillingList )
        np.save("fittingPara", fittingParaList )




    def factorFitting(self):
        """

        :return:
        """

        y =np.load("totalFluxNoiseHalfKDBSCAN.npy")
        x =self.rawBeamSize*self.smoothFactors

        #first dra this line

        fig = plt.figure(figsize=(10, 8))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        axFitting = fig.add_subplot(1,1, 1)

        pnDegree= 4
        z1 = np.polyfit(x, y, pnDegree )
        p1 = np.poly1d(z1)

        axFitting.plot( x,  p1(x) ,color='blue' ,lw=1 ,label="Polynomical fitting with degree {}".format(pnDegree ) )
        axFitting.scatter( x, y , s= 15, color='red',label="The noise RMS is 0.5 K" )

        ######
        axFitting.set_ylabel(r"Flux ($\rm K\ km\ s$$^{-1}$ $\mathrm{\Omega}_\mathrm{A}$)")
        axFitting.set_xlabel(r"Beam Size (arcmin)")
        print "Filling factor MWISP,", p1(self.rawBeamSize)/p1(0)
        print "Filling factor CfA,", p1(0.125*60)/p1(0)



        #draw CfA resolution
        axFitting.axvline(  x=8.5  , ls="--" ,  color='green' ,  label="CfA resolutuion (8.5 arcmin)" )
        axFitting.legend(loc= 3 )

        fig.tight_layout()
        plt.savefig("factorFitting.pdf", bbox_inches='tight')
        plt.savefig("factorFitting.png", bbox_inches='tight', dpi=600)


    def drawRMSwithBEAM(self):
        """
        Show show how rms change with beam size
        :return:
        """
        files=self.getSmFITSFileList()

        beamSize=self.rawBeamSize*self.smoothFactors

        if 0:
            rmsList=self.calRMSList( files  )
        else:
            rmsList=np.load("beamRMSQ1Local.npy")
        np.save("beamRMSQ1Local", rmsList )
        #draw

        fig = plt.figure(figsize=(10, 8))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        axRMS = fig.add_subplot(1,1, 1)



        ######
        axRMS.set_ylabel(r"RMS (K)")
        axRMS.set_xlabel(r"Beam Size (arcmin)")

        axRMS.plot(beamSize , rmsList   ,'o-',color='blue')
        #axRMS.plot(beamSize , rmsList   ,'o-',color='red', label="theoretical" )
        axRMS.plot(beamSize ,  rmsList[0]/self.smoothFactors   ,'o-',color='red', label="theoretical RMS" )

        #plot theretical


        #draw CfA resolution
        axRMS.axvline(  x=8.5  , ls="--" ,  color='green' ,  label="CfA resolutuion (8.5 arcmin)" )
        axRMS.legend(loc= 1 )

        fig.tight_layout()
        plt.savefig("rmsWithBeam.pdf", bbox_inches='tight')
        plt.savefig("rmsWithBeam.png", bbox_inches='tight', dpi=600)


    def pipeLineTestabsK(self,downToK=1.28, recalTotalFlux= False):
        """
        used to compare different method to calculate flux, in this pipeline, the noise is fixed to be 0.5 K, then change the beam size, what we would see?
        :return:
        """
        #first


        if recalTotalFlux:

            TBList,tbFileList= self.getAbsNoiseCleanTBLIst()
            totalFluxDBSCAN=self.getTotalFlux(TBList)
            np.save("totalFluxNoiseHalfKDBSCAN" , totalFluxDBSCAN)



            fitsFileList = self.getAbsNoiseFITSName()

            fluxListDirectSum = self.sumVoxels(  fitsFileList,  downToK=downToK  )

            print "Saving total flux, "
            np.save("totalFluxNoiseHalfKSum" , fluxListDirectSum)

            for eachF in fitsFileList:

                continue # Do not do this a second time

                self.cleanFITSsigma2(eachF,removeFITS=False)



        else:

            totalFluxDBSCAN=np.load("totalFluxNoiseHalfKDBSCAN.npy")
            fluxListDirectSum=np.load("totalFluxNoiseHalfKSum.npy")


        if 0:

            eqSigmaList=[]

            for i in range(len(fitsFileList)):

                fitsName= fitsFileList[i]
                tbFileName=  tbFileList[i]
                eqDownToSigma=self.getEqualCutOff(fitsName, tbFileName)/self.MWISPrmsCO12

                eqSigmaList.append(eqDownToSigma)



            print eqSigmaList
            print "The mean sigma is,", np.mean(eqSigmaList)


        #total flux ,simply sum all voxels

        #plot those to flux


        fig = plt.figure(figsize=(10, 8))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]


        axBeam = fig.add_subplot(1,1, 1)
        axBeam.plot( self.smoothFactors*self.rawBeamSize, totalFluxDBSCAN,     'o-',   color=  'blue'  ,markersize=2 ,lw=1,label= "FLux by DBSCAN clouds"  )


        downToSimg=downToK/self.MWISPrmsCO12
        axBeam.plot( self.smoothFactors*self.rawBeamSize, fluxListDirectSum,     'o-',   color=  'green'  ,markersize=2 ,lw=1,label= r"FLux by summing voxels above {}$\sigma$".format(downToSimg)  )

        at = AnchoredText(r"The noise rms is {} K (1$\sigma$)".format(self.MWISPrmsCO12), loc=1, frameon=False)
        axBeam.add_artist(at)



        #draw with color

        #self.drawBeamEffectSingle(axBeam, noiseFactor=1.0 )
        #self.drawBeamEffectSingle(axBeam, noiseFactor=1.5  )


        #self.drawNoiseEffectSingle(axSens,1.5)

        axBeam.set_ylabel(r"Flux ($\rm K\ km\ s$$^{-1}$ $\mathrm{\Omega}_\mathrm{A}$)")
        axBeam.set_xlabel(r"Beam Size (arcmin)")

        #self.drawColorBarNoiseFactor(axBeam)
        #draw cfa beam factor

        axBeam.axvline(  x=8.5  ,  color='purple' ,  label="CfA resolutuion (8.5 arcmin)" )
        axBeam.legend(loc= 3 )

        fig.tight_layout()
        plt.savefig("effectBeamFixNoise.pdf", bbox_inches='tight')
        plt.savefig("effectBeamFixNoise.png", bbox_inches='tight', dpi=600)





    def getFITSBySF(self, SF=1 , calCode=None):

        if calCode is None:
            calCode=self.calCode


        SF=np.float(SF)
        searchStr="{}*SmFactor_{}.fits".format(self.calCode,SF)
        searchStr= os.path.join( self.tmpPath,searchStr )
        foundFiles= glob.glob(searchStr)

        if len( foundFiles )==1:
            return foundFiles[0]

        else:
            print searchStr
            print "Please check you data!"
            sys.exit()



    def drawNoiseEffectSingle(self,ax,smFactor):

        tbList,noiseArray=self.getTBFixSmFactor(smFactor=smFactor)

        totalFlux=self.getTotalFlux(tbList)


        #better convert noise factor to Kelvin
        #get the noise factor

        fitsRawFile=self.getFITSBySF(SF=smFactor   )

        factor1RMS=self.calRMSSingle(fitsRawFile, onlyValue=True)
        print "The rms for factor 1 is, ", factor1RMS

        plotColor = self.smoothColor.to_rgba(smFactor)

        ax.plot(noiseArray*factor1RMS, totalFlux  ,'o-',   color=  plotColor  ,markersize=2 ,lw=1)




    def drawfluxChange(self, cutoff=0.5):

        fileList,smoothFactors=self.getSmFITSFileList()

        fluxList_1 =[]
        fluxList_2 =[]
        fluxList_3 =[]
        fluxList_4 =[]
        fluxList_5 =[]
        fluxList_6 =[]
        fluxList_7 =[]

        omega = np.deg2rad(0.5/60)*np.deg2rad(0.5/60)/4./np.log(2)

        rms=0.5

        for eachFITS in fileList:
            print "Processing ",eachFITS
            data,head=doFITS.readFITS(eachFITS)

            goodValues= data>=1*rms
            data[~goodValues]=0
            fluxList_1.append( np.sum(data)*self.getVelResolution()*omega  )

            data[data<2*rms]=0
            fluxList_2.append( np.sum(data)*self.getVelResolution()*omega  )

            data[data<3*rms]=0
            fluxList_3.append( np.sum(data)*self.getVelResolution()*omega  )


            data[data<4*rms]=0
            fluxList_4.append( np.sum(data)*self.getVelResolution()*omega  )


            data[data<5*rms]=0
            fluxList_5.append( np.sum(data)*self.getVelResolution()*omega  )

            data[data<6*rms]=0
            fluxList_6.append( np.sum(data)*self.getVelResolution()*omega  )

            data[data<7*rms]=0
            fluxList_7.append( np.sum(data)*self.getVelResolution()*omega  )






        #plot figures

        fig = plt.figure(figsize=(8, 6))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 13, 'serif': ['Helvetica']})
        axTotalFlux = fig.add_subplot(1, 1, 1)

        axTotalFlux.plot(smoothFactors,fluxList_1,'o-',lw=1,label=r"1$\sigma$ (0.5 K)")
        axTotalFlux.plot(smoothFactors,fluxList_2,'o-',lw=1,label=r"2$\sigma$ (1.0 K)")
        axTotalFlux.plot(smoothFactors,fluxList_3,'o-',lw=1,label=r"3$\sigma$ (1.5 K)")
        axTotalFlux.plot(smoothFactors,fluxList_4,'o-',lw=1,label=r"4$\sigma$ (2.0 K)")

        axTotalFlux.plot(smoothFactors,fluxList_5,'o-',lw=1,label=r"5$\sigma$ (2.5 K)")
        axTotalFlux.plot(smoothFactors,fluxList_6,'o-',lw=1,label=r"6$\sigma$ (3.0 K)")
        axTotalFlux.plot(smoothFactors,fluxList_7,'o-',lw=1,label=r"7$\sigma$ (3.5 K)")


        axTotalFlux.legend(loc=1)

        axTotalFlux.set_ylabel(r"Total flux ($\rm K\ km\ s$$^{-1}$ $\Omega_\mathrm{A}$)")
        axTotalFlux.set_xlabel(r"Smooth factor")


        plt.savefig( "G2650smoothFactorFlux.png", bbox_inches='tight', dpi=300)
        plt.savefig(  "G2650smoothFactorFlux.pdf", bbox_inches='tight')


    def addNoiseToWMISPNoise(self  ):
        """

        :return:
        """

        smoothFiles= self.getSmFITSFileList()

        for eachFile in  smoothFiles:
            self.addNoiseSingle(  eachFile,  absoluteNoise=self.MWISPrmsCO12 )



    def getRawBeamTBByCalcode(self,calCode,smFactor=1):

        searchStr= "{}*SmFactor_{}*absKdbscanS2P4Con1_Clean.fit".format(calCode, float(smFactor) )

        testSearch  = glob.glob( self.tmpPath+searchStr )


        return   testSearch[0]

    def getCleanFITSName(self,calCode,smFactor):

        searchStr= "{}*SmFactor_{:.1f}*absKdbscanS2P4Con1_Clean.fits".format(calCode,float(smFactor))
        testSearch  = glob.glob( self.tmpPath+searchStr )


        return   testSearch[0]

    def getRawCOFITS(self,calCode=None):
        """

        :param calCode:
        :return:
        """
        if calCode!=None:

            return self.dataPath+calCode+".fits"
        else:
            return self.dataPath+self.calCode+".fits"




    def getFluxListByIDTest2(self,calCode,ID):
        """

        :param ID:
        :return:
        """

        fluxList=[]

        rawCOFITS=self.getRawCOFITS(calCode)

        CODataRaw, COHeadRaw = doFITS.readFITS( rawCOFITS )

        cleanFITSRawBeam = self.getCleanFITSName(calCode, 1)

        #step 1, get the
        data,head=doFITS.readFITS(cleanFITSRawBeam)
        cloudIndex= np.where( data==ID )

        rawCOValues=  CODataRaw[cloudIndex]


        for eachSm in self.smoothFactors:

            print eachSm
            cleanFITSSm =  self.getCleanFITSName(calCode, eachSm )

            dataSm, headSm = doFITS.readFITS( cleanFITSSm )

            #smLabels=  dataSm[cloudIndex]

            indexSelect = np.logical_and(  data==ID , dataSm>0 )

            coValues= CODataRaw[indexSelect]    #rawCOValues[smLabels>0]
            flux=np.sum(coValues)*self.getVelResolution()
            print flux
            fluxList.append(flux)



        return fluxList


    def getCleanDataList(self,calCode):

        dataArray=[]

        for eachSm in self.smoothFactors:
            
            cleanFITSSm =  self.getCleanFITSName(calCode, eachSm )

            dataSm, headSm = doFITS.readFITS( cleanFITSSm )
            dataArray.append( dataSm )

        return dataArray


    def getFluxListByIDSigmaCut(self,  labelSets, smCODataList  ,ID ,sigmaCut=2.6 ):
        """

        :param ID:
        :return:
        """

        fluxList=[]

        #rawCOFITS=self.getRawCOFITS(calCode)

        #CODataRaw, COHeadRaw = doFITS.readFITS( rawCOFITS )

        #step 1, get the
        #data,head=doFITS.readFITS(cleanFITSRawBeam)
        cloudIndex=  self.getIndices(labelSets,ID)  #np.where( cleanDataSM1==ID )

        #rawCOValues=  CODataRaw[cloudIndex]

        noise=self.guessNoiseByCode()

        for  dataSm in smCODataList:

            coValues=   dataSm[cloudIndex]
            coValues= coValues[coValues>=sigmaCut* noise ]
            flux=np.sum(coValues)*self.getVelResolution()
            fluxList.append(flux)


        del cloudIndex
        del coValues


        return fluxList


    def addCutFFColnames(self,cleanTB):
        ffTB=cleanTB.copy()
        ffTB[self.cutFFffMWISPCol] = cleanTB["peak"] * 0
        ffTB[self.cutFFffMWISPCol].unit=""

        
        ffTB[self.cutFFffCfACol] = cleanTB["peak"] * 0
        ffTB[self.cutFFffCfACol].unit=""

        ffTB[self.cutFFaCol] = cleanTB["peak"] * 0
        ffTB[self.cutFFaCol].unit=""

        ffTB[self.cutFFbCol] = cleanTB["peak"] * 0
        ffTB[self.cutFFcCol] = cleanTB["peak"] * 0
        ffTB[self.cutFFaErrorCol] = cleanTB["peak"] * 0
        ffTB[self.cutFFbErrorCol] = cleanTB["peak"] * 0
        ffTB[self.cutFFcErrorCol] = cleanTB["peak"] * 0

        ffTB[self.cutFFbCol].unit=""
        ffTB[self.cutFFcCol].unit=""
        ffTB[self.cutFFaErrorCol].unit=""
        ffTB[self.cutFFbErrorCol].unit=""
        ffTB[self.cutFFcErrorCol].unit=""


        return ffTB

    def addFFColnames(self,cleanTB):
        ffTB=cleanTB.copy()
        ffTB[self.ffMWISPCol] = cleanTB["peak"] * 0
        ffTB[self.ffCfACol] = cleanTB["peak"] * 0
        ffTB[self.aCol] = cleanTB["peak"] * 0
        ffTB[self.bCol] = cleanTB["peak"] * 0
        ffTB[self.cCol] = cleanTB["peak"] * 0

        ffTB[self.aErrorCol] = cleanTB["peak"] * 0
        ffTB[self.bErrorCol] = cleanTB["peak"] * 0
        ffTB[self.cErrorCol] = cleanTB["peak"] * 0

        return ffTB

    def addCutOffFluxColnames(self, cleanTB ):
        """
        adding columns to record cutoff flux
        :param cleanTB:
        :return:
        """
        for eachCut in self.cutoffFactors:
 
            colName,colPixName=self.getFluxColNameCutoff(eachCut)

            cleanTB[colName]=cleanTB["peak"]*0
            cleanTB[colPixName]=cleanTB["peak"]*0
            cleanTB[colPixName].unit=""





    def getIndices(self,labelSets, choseID ,return2D=False):
        Z0, Y0, X0, values1D =labelSets
        cloudIndices = np.where(values1D == choseID)

        cX0 = X0[cloudIndices]
        cY0 = Y0[cloudIndices]
        cZ0 = Z0[cloudIndices]

        if not return2D:

            return tuple([cZ0, cY0, cX0])
        else:
            return tuple([cZ0, cY0, cX0]), tuple([  cY0, cX0])

            
    def getCloudCube(self,ID):

        searCubePath=self.checkCloudCubeSavePath()

        searchStr=os.path.join(searCubePath,self.calCode+"*{}*.fits".format(ID) )



        searchResults = glob.glob( searchStr )

        return searchResults[0]





    def drawInMap(self,calCode,ID):
        """

        get the sub cube and draw the integration map

        :param calCode:
        :param ID:
        :return:
        """

        self.calCode = calCode

        cloudCube=self.getCloudCube(ID)


        dataCO,head=doFITS.readFITS(cloudCube)

        intMap=np.sum(dataCO,axis=0)

        ##########################
        #draw integration map

        fig = plt.figure(1, figsize=(15,8))
        rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica'], "size": 15 })
        # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

        rc('text', usetex=True)
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]
        WCSCrop=WCS(head)
        axCO = pywcsgrid2.subplot(111, header=WCSCrop)


        cmapCO = plt.cm.bone
        cmapCO.set_bad('black')
        axCO.imshow(np.sqrt(intMap), origin='lower', cmap=cmapCO, vmin=0, vmax=3, interpolation='none')

        axCO.set_ticklabel_type("absdeg", "absdeg")
        axCO.axis[:].major_ticks.set_color("w")
        axCO.set_xlabel(r"Galactic Longitude ($^{\circ}$)")
        axCO.set_ylabel(r"Galactic Latitude ($^{\circ}$)")
        saveName=self.intFigurePath+"{}_{}CloudInt".format(calCode,ID)
        plt.savefig("compareSCIMESPara.pdf" , bbox_inches="tight")
        plt.savefig(saveName+".png" , bbox_inches="tight", dpi=600)





    def calFFByID(self,calCode,ID, cleanTB, drawFigure=True, useSigmaCut=True  ):

        """
        calculate filling factor, for only one  specific cloud
        ususlly this is used to check filling factor for a specific cloud
        :param calCode:
        :param drawFigure:
        :param useSigmaCut:
        :return:
        """

        self.calCode= calCode

        print calCode
        #TBName =  self.getSmoothAndNoiseFITSSingle(getCleanTBFile=True)#self.getRawBeamTBByCalcode(calCode)
        #cleanTB = Table.read(TBName)

        #ffTB=self.addFFColnames( cleanTB )

        rowID= self.getRowByID(cleanTB,ID)

        self.calFFByRow(rowID,drawFigure=True)




    def getRowByID(self,TB,ID):
        """

        :param TB:
        :param id:
        :return:
        """
        #fined the index of ID

        idCol=TB["_idx"]

        indexID = np.where(idCol==ID)
        indexID=indexID[0][0]

        return TB[indexID]




    def getFluxColName(self,smFactor):
        """
        The pix name is is to record the tao pix Number of the flux, which is used to estimate the error of the total flux
        :param smFactor:
        :return:
        """
        return "fluxSM{:.1f}".format( smFactor ), "pixSM{:.1f}".format( smFactor )



    def getFluxColNameCutoff(self,cutoffFactor):
        """
        The pix name is is to record the tao pix Number of the flux, which is used to estimate the error of the total flux
        :param smFactor:
        :return:
        """
        return "fluxCut{:.1f}".format( cutoffFactor ), "pixCut{:.1f}".format( cutoffFactor )




    #def getFluxColNameNoiseChange(self,cutoffFactor):
        #"""
        #The pix name is is to record the tao pix Number of the flux, which is used to estimate the error of the total flux
        #:param smFactor:
        #:return:
        #"""
        #return "fluxNC{:.1f}".format( cutoffFactor ), "pixNC{:.1f}".format( cutoffFactor )








    def getSmoothFluxCol(self,smFITS, TB,  labelSets, sigmaCut=2.6 ):
        """
        :return:
        """

        ####
        print "Extracting flux from ", smFITS
        dataSm,headSm= doFITS.readFITS(smFITS)

        smFactor = self.getSmoothFactor(smFITS)
        colName,colPixName=self.getFluxColName(smFactor)

        TB[colName]=TB["peak"]*0
        TB[colPixName]=TB["peak"]*0

        widgets = ['Extracting flux:', Percentage(), ' ', Bar(marker='0', left='[', right=']'), ' ', ETA(), ' ',
                   FileTransferSpeed()]  # see docs for other options
        pbar = ProgressBar(widgets=widgets, maxval=len(TB))
        pbar.start()


        i=0
        noise=self.getMeanRMS()
        for eachRow in TB:
            i=i+1
            #gc.collect()

            ID = eachRow["_idx"]

            cloudIndex = self.getIndices(labelSets, ID)  # np.where( cleanDataSM1==ID )
            coValues=   dataSm[cloudIndex]
            coValues= coValues[coValues>=sigmaCut* noise ]
            fluxID=np.sum(coValues)*self.getVelResolution()

            eachRow[colName] = fluxID
            eachRow[colPixName] =  len(coValues)

            pbar.update(i)

        pbar.finish()


    def getVelResolution(self,calCode=None):

        """
        return the velocity resolution in units of km/s
        :param calCode:
        :return:
        """

        if calCode==None:
            testCode=self.calCode
        else:
            testCode = calCode

        if testCode=="":
            print "Please set you calculation code!"
            aaaaaaaaaaaaa


        isRawCode = self.isRaw(testCode)

        if isRawCode:
            if "12" in testCode:
                return 0.158737644553

            if "13" in testCode:
                return 0.166040420532

            if "18" in testCode:
                return 0.166674420237

        else: #if not raw, must be 0.2 km/s
            return 0.2









    def getFluxList(self,row, useCutoff=False ):
        """

        :param row:
        :return:
        """


        if useCutoff:

            fluxList = []
            fluxErrorList = []

            meanRMS = self.getMeanRMS()
            
            
            
            for eachCutoff in self.cutoffFactors :
                colName, colPixName = self.getFluxColNameCutoff(eachCutoff)
                fluxList.append(row[colName])  # the flux is already K km/s
                totalVox = row[colPixName]

                totalVox = max([1, totalVox])  # assign 1 pixel error to 0 flux

                eRRor = np.sqrt(totalVox) * meanRMS * self.getVelResolution()

                fluxErrorList.append(eRRor)

            return np.asarray(fluxList), fluxErrorList




        fluxList=[]
        fluxErrorList=[]

        meanRMS=self.getMeanRMS()
        for eachSm in self.smoothFactors:
            colName,colPixName=self.getFluxColName(eachSm)
            fluxList.append( row[colName]) #the flux is already K km/s
            totalVox=row[colPixName]

            totalVox=max([1,totalVox]) #  assign 1 pixel error to 0 flux

            eRRor= np.sqrt(totalVox)*meanRMS*self.getVelResolution()



            fluxErrorList.append( eRRor)



        return np.asarray(fluxList), fluxErrorList






    def calFFByRow(self,eachRow,drawFigure=False):
        """
        calculate filling factors by row
        :param tbRow:
        :return:
        """
        processBeam = self.getBeamSize() #in arcmin

        beamArray= self.smoothFactors* processBeam

        fluxList = self.getFluxList(eachRow)
        ID = eachRow[self.idCol]
        wmsipFilling, cfaFilling, fittingParaAndError = self.getFillingFactorAndDraw(beamArray, fluxList, calID=ID,   drawFigure=drawFigure)
        return wmsipFilling


    def calculateFillingFactorCutoff(self,TB,drawFigure=False,inputID=None):

        """
        :param TB:
        :return:
        """
        #TB=Table.read(TBFile)
        #saveName="cutFF_"+TBFile

        TB=self.addCutFFColnames(TB)

        cutoffArray= self.cutoffFactors* self.getMeanRMS()

        #add progress bar
        widgets = ['Calculating cutoff filling factors: ', Percentage(), ' ', Bar(marker='0', left='[', right=']'), ' ', ETA(), ' ',
                   FileTransferSpeed()]  # see docs for other options
        pbar = ProgressBar(widgets=widgets, maxval=len(TB))
        pbar.start()

        i=0
        for eachRow in TB:
            i=i+1
            pbar.update(i)
            fluxList,fluxError = self.getFluxList(eachRow,useCutoff=True)


            ID=eachRow[ self.idCol ]
            if inputID!=None and ID==inputID: #for debug

                print fluxList
                wmsipFilling, cfaFilling, fittingParaAndError = self.getFillingFactorAndDraw(cutoffArray, fluxList,   fluxError, calID=ID,   drawFigure=True)
                print "filling factor", wmsipFilling
                print "Parameters", fittingParaAndError[0]
                print "Para error", fittingParaAndError[1]

                return

            #crop


            wmsipFilling, cfaFilling,  fittingParaAndError=self.getFillingFactorAndDraw(cutoffArray,fluxList,fluxError,calID=ID,drawFigure=drawFigure)

            para, paraError = fittingParaAndError
            eachRow[self.cutFFffMWISPCol]=wmsipFilling
            eachRow[self.cutFFffCfACol]=cfaFilling

            eachRow[self.cutFFaCol]=para[0]
            eachRow[self.cutFFbCol]=para[1]
            eachRow[self.cutFFcCol]=para[2]


            eachRow[self.cutFFaErrorCol]=paraError[0]
            eachRow[self.cutFFbErrorCol]=paraError[1]
            eachRow[self.cutFFcErrorCol]=paraError[2]






        pbar.finish()
        #TB.write( saveName , overwrite=True )
        return TB







    def calculateFillingFactor(self,TBFile,drawFigure=False,inputID=None):

        """
        :param TB:
        :return:
        """
        TB=Table.read(TBFile)
        saveName="fillingFactor_"+TBFile
        processBeam = self.getBeamSize() #in arcmin

        beamArray= self.smoothFactors* processBeam

        #add progress bar
        widgets = ['Calculating filling factors: ', Percentage(), ' ', Bar(marker='0', left='[', right=']'), ' ', ETA(), ' ',
                   FileTransferSpeed()]  # see docs for other options
        pbar = ProgressBar(widgets=widgets, maxval=len(TB))
        pbar.start()

        i=0
        for eachRow in TB:
            i=i+1
            pbar.update(i)
            fluxList,fluxError = self.getFluxList(eachRow)



            ID=eachRow[ self.idCol ]
            if inputID!=None and ID==inputID: #for debug

                print fluxList
                wmsipFilling, cfaFilling, fittingParaAndError = self.getFillingFactorAndDraw(beamArray, fluxList,   fluxError, calID=ID,   drawFigure=True)
                print "filling factor", wmsipFilling
                print "Parameters", fittingParaAndError[0]
                print "Para error", fittingParaAndError[1]

                return

            wmsipFilling, cfaFilling,  fittingParaAndError=self.getFillingFactorAndDraw(beamArray,fluxList,fluxError,calID=ID,drawFigure=drawFigure)



            para, paraError = fittingParaAndError
            eachRow[self.ffMWISPCol]=wmsipFilling
            eachRow[self.ffCfACol]=cfaFilling

            eachRow[self.aCol]=para[0]
            eachRow[self.bCol]=para[1]
            eachRow[self.cCol]=para[2]


            eachRow[self.aErrorCol]=paraError[0]
            eachRow[self.bErrorCol]=paraError[1]
            eachRow[self.cErrorCol]=paraError[2]
        pbar.finish()
        TB.write( saveName , overwrite=True )
        return saveName


    def removeScutumBadChannels(self,TB):
        """
        The TB should be cloud catalog of the Scutum Arm
        :param TB:
        :return:
        """

        #scutum arm, at about 96 km/s, for CO13 an CO18, area bad channels

        #rawOutCO12:
        #bad v channel, -24.286

        #ppv1 v[-25, -23],  box(47.4879697,2.9315877,15379.103 ,14716.975 ,0)
        #ppv2 v[-25, -23], box(45.2827582,-1.6997356,8640.000 ,8460.000 ,0)


        #rawOutCO13
        #ppv1 v, less than -69 km/s, remove all

        # rawOutCO18
        # ppv1 v, less than -69 km/s, remove all




    def getFluxListForEachCloudNoiseChange(self,calCode=None, drawFigure=False, useSigmaCut=True, calAllCloud=True ):

        """
        Due to the memory problem, this part of code need to be revise

        #step 1 produce flux table,
        #step 2, calculate fff with the flux table
        :param calCode:
        :param drawFigure:
        :param useSigmaCut:
        :param calAllCloud:
        :return:
        """
        if calCode!=None:
            self.calCode= calCode



        TBName = self.getSmoothAndNoiseFITSSingle(smFactor=1.0,noiseFactor=0.0,getCleanTBFile=True)


        cleanTB = Table.read(TBName)
        ffTB=self.addFFColnames( cleanTB )

        #rawCOFITS=self.getRawCOFITS(calCode)
        #CODataRaw, COHeadRaw = doFITS.readFITS( rawCOFITS )

        cleanFITSRawBeam = self.getSmoothAndNoiseFITSSingle(smFactor=1.0,noiseFactor=0.0,  getCleanFITS =True)
        cleanDataSM1,head=doFITS.readFITS(cleanFITSRawBeam)

        clusterIndex1D = np.where(cleanDataSM1 > 0)
        clusterValue1D = cleanDataSM1[clusterIndex1D]

        Z0, Y0, X0 = clusterIndex1D
        labelSets=[Z0, Y0, X0, clusterValue1D ]



        #the next step is to extract flux
        #getSmoothFluxColNoiseChange
        #allSmoothFiles = self.getSmoothListFixNoise(noiseFactor=0.)
        allNoiseFiles= self.getSmoothListFixSmooth(smoothFactor=1.0)
        for eachSmFile in allNoiseFiles:

            self.getSmoothFluxColNoiseChange(eachSmFile,ffTB,labelSets )


        ffTB.write("fluxTBNoiseChange_" + os.path.basename(TBName) ,overwrite=True )

        #step, fitting filling factor


        return
        ####

    def getFluxListForEachCloud(self,calCode=None, drawFigure=False, useSigmaCut=True, calAllCloud=True ):

        """
        Due to the memory problem, this part of code need to be revise

        #step 1 produce flux table,
        #step 2, calculate fff with the flux table
        :param calCode:
        :param drawFigure:
        :param useSigmaCut:
        :param calAllCloud:
        :return:
        """
        if calCode!=None:
            self.calCode= calCode



        TBName = self.getSmoothAndNoiseFITSSingle(smFactor=1.0,noiseFactor=0.0,getCleanTBFile=True)



        cleanTB = Table.read(TBName)

        ffTB=self.addFFColnames( cleanTB )

        #rawCOFITS=self.getRawCOFITS(calCode)
        #CODataRaw, COHeadRaw = doFITS.readFITS( rawCOFITS )

        cleanFITSRawBeam = self.getSmoothAndNoiseFITSSingle(smFactor=1.0,noiseFactor=0.0,  getCleanFITS =True)
        cleanDataSM1,head=doFITS.readFITS(cleanFITSRawBeam)

        clusterIndex1D = np.where(cleanDataSM1 > 0)
        clusterValue1D = cleanDataSM1[clusterIndex1D]

        Z0, Y0, X0 = clusterIndex1D
        labelSets=[Z0, Y0, X0, clusterValue1D ]



        #the next step is to extract flux

        allSmoothFiles = self.getSmoothListFixNoise(noiseFactor=0.)

        for eachSmFile in allSmoothFiles:

            self.getSmoothFluxCol(eachSmFile,ffTB,labelSets )


        ffTB.write("fluxTB_" + os.path.basename(TBName) ,overwrite=True )

        #step, fitting filling factor


        return
        ####



        #remove small size clouds
        if not calAllCloud:
            area=ffTB["area_exact"]

            size= 2*np.sqrt(area/np.pi)

            ffTB=ffTB[size>= self.smoothFactors[-1]*self.rawBeamSize]



        #CODataRaw, cleanDataSM1, cleanDataList,



        if useSigmaCut:
            cleanDataList=self.getNoiseDataList()

        else:
            cleanDataList=self.getCleanDataList(self.calCode)

        widgets = ['Calculating filling factors: ', Percentage(), ' ', Bar(marker='0', left='[', right=']'), ' ', ETA(), ' ',
                   FileTransferSpeed()]  # see docs for other options
        pbar = ProgressBar(widgets=widgets, maxval=len(ffTB))
        pbar.start()


        i=0

        clusterIndex1D = np.where(cleanDataSM1 > 0)
        clusterValue1D = cleanDataSM1[clusterIndex1D]

        Z0, Y0, X0 = clusterIndex1D
        labelSets=[Z0, Y0, X0, clusterValue1D ]


        for eachR in ffTB:
            i=i+1
            #gc.collect()

            ID=eachR["_idx"]
            pbar.update(i)

            self.getFillingFactorByCloudID(CODataRaw, labelSets,cleanDataList,self.calCode, ID, saveRow=eachR, drawFigure=drawFigure, useSigmaCut=useSigmaCut)


        pbar.finish()
            #break
        if calAllCloud:
            ffTB.write( self.calCode+"FillingFactorTBAll.fit",overwrite=True )

        else:
            ffTB.write( self.calCode+"FillingFactorTB.fit",overwrite=True )


    def testFunctionForm(self):
        """

        :return:
        """

        fig = plt.figure(figsize=(10, 8))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 17, 'serif': ['Helvetica']})


        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        axPlot = fig.add_subplot(1,1, 1)

        X=np.arange(0.1,10,0.001)

        Y = ffAndSizeFunc(X,0.5)

        axPlot.plot(X,Y,'b-',lw=1)

        plt.savefig(  "functionTest.png" , bbox_inches='tight', dpi=300)


    def cloudStat(self, drawTBFile=None):
        """
        draw statistics of Angular area, Flux, Vrms, and Peak
        :return:
        """

        #four figures

        drawTBFile= self.getSmoothAndNoiseFITSSingle(getCleanTBFile=True)

        drawTB=Table.read( drawTBFile )

        print "Total number molecular clouds, ", len(drawTB)

        fig = plt.figure(figsize=(20, 16))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 17, 'serif': ['Helvetica']})


        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]



        axArea=fig.add_subplot(2,2, 1)

        #draw hist of ara

        self.drawAreaHist(axArea,drawTB)
        #mark tag here

        at = AnchoredText(self.calCode , loc=1, frameon=False)
        axArea.add_artist(at)



        axFlux=fig.add_subplot(2,2, 2)

        self.drawFluxHist(axFlux,drawTB)

        axPeak=fig.add_subplot(2,2, 3)
        self.drawPeakHist(axPeak,drawTB)
        axVel=fig.add_subplot(2,2, 4)
        self.drawVrmsHist(axVel,drawTB)

        plt.savefig( self.statPath+"{}CloudStatistic.png".format(self.calCode  ), bbox_inches='tight', dpi=600)


    def drawAreaHist(self,ax,drawTB):

        """
        draw histgram of area
        :param ax:
        :param TB:
        :return:
        """
        ax.set_yscale('log')
        ax.set_xscale('log')
        
        
        areaEdges = np.linspace(0.25 / 3600., 150, 10000)
        areaCenter = self.getEdgeCenter(areaEdges)
        #first, get cloud catalog

        binN, binEdges = np.histogram(drawTB["area_exact"] / 3600., bins=areaEdges)

        ax.plot(areaCenter[binN > 0], binN[binN > 0], color= "blue", linestyle='-', markersize=3,
                lw=1, marker='o', markeredgewidth=0.6, fillstyle='full', alpha=0.8, )
        ax.set_ylabel(r"Number of molecular clouds")


        ax.set_xlabel(r"Angular area (deg$^2$)")


    def drawFluxHist(self,ax,drawTB):
        """
        draw hist of flux
        :param ax:
        :param drawTB:
        :return:
        """
        ax.set_yscale('log')
        ax.set_xscale('log')


        fluxEdges = np.linspace(8, 1e5, 1000)
        fluxCenter = self.getEdgeCenter( fluxEdges )

        sum = drawTB["sum"] * self.getVelResolution()  # K km/s
        binN, binEdges = np.histogram(sum, bins=fluxEdges)

        ax.plot(fluxCenter[binN > 0], binN[binN > 0], color="green", linestyle='-',
                markersize=3, lw=1, marker='o', markeredgewidth=0.5, fillstyle='full', alpha=0.8, )

        strXlabel = r"Flux ($\rm K\ km\ s$$^{-1}$ $\mathrm{\Omega}_{\rm A}$)"

        ax.set_xlabel( strXlabel )

        ax.set_ylabel(r"Number of molecular clouds")


    def drawPeakHist(self,ax,drawTB):
        """
        draw histgram of peak distributuion
        :return:
        """

        velEdges = np.linspace(0, 40, 400)
        velCenter = self.getEdgeCenter(velEdges)

        vDisperse = drawTB["peak"]

        binN , binEdges  = np.histogram(vDisperse, bins=velEdges)

        stepa = ax.step(velCenter, binN, lw= 1 ,  linestyle='-', color= "black"   )

        ax.set_ylabel(r"Number of clusters")
        ax.set_xlabel(r"Peak brightness temperature (K)" )

        ax.set_xlim(0, 10)





    def drawVrmsHist(self,ax,drawTB):

        """
        #draw histgram of Vrms
        :param ax:
        :param drawTB:
        :return:
        """


        velEdges = np.linspace(0, 15, 300)
        velCenter = self.getEdgeCenter(velEdges)
        #################################

        vDisperse = drawTB["v_rms"]

        binN , binEdges  = np.histogram(vDisperse, bins=velEdges)



        ax.step(velCenter, binN, lw= 1 ,  linestyle='-', color= 'purple'   )
        ax.set_ylabel(r"Number of molecular clouds")
        ax.set_xlabel(r"Velocity dispersion ($\rm km\ s$$^{-1}$)")
        ax.set_xlim([0, 3])

    def getEdgeCenter(self, edges):

        areaCenters = (edges[1:] + edges[0:-1]) / 2.

        return areaCenters

    def getFillingErrorAndError(self, beamSize , ffTB, useCutoff=False ):

        """

        :param ffTB:
        :return:
        """
        print ffTB.colnames

        if useCutoff:
            x = beamSize  # arcmin
            a = ffTB[self.cutFFaCol]
            b = ffTB[self.cutFFbCol]
            c = ffTB[self.cutFFcCol]

            stdA = ffTB[self.cutFFaErrorCol]
            stdB = ffTB[self.cutFFbErrorCol]
            stdC = ffTB[self.cutFFcErrorCol]

            fx = ffFunction(x, a, b, c)

            f0 = ffFunction(0, a, b, c)
            varX = stdA ** 2 * (dya(x, a, b, c)) ** 2 + stdB ** 2 * (dyb(x, a, b, c)) ** 2 + stdC ** 2

            var0 = stdA ** 2 * (dya(0, a, b, c)) ** 2 + stdB ** 2 * (dyb(0, a, b, c)) ** 2 + stdC ** 2

            varFF = varX / f0 ** 2 + (-fx / f0 ** 2) ** 2 * var0

            return fx / f0, np.sqrt(varFF)


        x=  beamSize #arcmin
        a=ffTB[self.aCol]
        b=ffTB[self.bCol]
        c=ffTB[self.cCol]

        stdA = ffTB[ self.aErrorCol ]
        stdB = ffTB[ self.bErrorCol ]
        stdC = ffTB[ self.cErrorCol ]

        fx=ffFunction(x,a,b,c)

        f0= ffFunction(0,a,b,c)
        varX = stdA**2* ( dya(x,a,b,c) )**2 + stdB**2* ( dyb(x,a,b,c) )**2 + stdC**2


        var0=  stdA**2* ( dya(0,a,b,c) )**2 + stdB**2* ( dyb(0,a,b,c) )**2 + stdC**2


        varFF=varX/f0**2 + (-fx/f0**2)**2*var0

        return  fx/f0  ,  np.sqrt( varFF )
        #beamSize


    def addMWISPFFerror(self, TB  ):
        """

        :return:
        """

        #TB = Table.read( TBName )
        #print len(TB)

        #TB=TB[TB[self.ffMWISPCol]>0]
        #print len(TB)

        mwispFF,mwispFFError= self.getFillingErrorAndError(self.getBeamSize(), TB   )


        #relativeError= mwispFFError/mwispFF


        #goodFF =  TB[ relativeError <0.3]

        #print mwispFFError

        TB[self.ffMWISPErrorCol] = mwispFFError

        return TB


    def addMWISPFFerrorCutoff(self, TB  ):
        """

        :return:
        """

        mwispFF,mwispFFError= self.getFillingErrorAndError(self.getBeamSize(), TB  ,useCutoff=True )


        TB[self.cutFFffMWISPErrorCol] = mwispFFError

        return TB


    def addCfaFFerror(self, TB ):
        """

        :return:
        """


        cfaFF, cfaFFError= self.getFillingErrorAndError(8.5, TB   )


        TB[self.ffCfAErrorCol] = cfaFFError

        return TB




    def getCloudWithGoodFF(self,ffTB):

        ffTB=ffTB[ ffTB[self.ffMWISPCol]>0 ]

        errorA = ffTB[self.aErrorCol]/ffTB[self.aCol]
        errorB = ffTB[self.bErrorCol]/ffTB[self.bCol]
        errorC = ffTB[self.cErrorCol]/ffTB[self.cCol]

        errorThreshold = 0.3
        selectCriteria=np.logical_and(errorA< errorThreshold ,errorB< errorThreshold )

        selectCriteria=np.logical_and( selectCriteria  , errorC< errorThreshold )



        ffTB=ffTB[ selectCriteria ]


        return  ffTB


    def getEdgeClouds(self, TB):
        """

        :param TB:
        :return:
        """

        selectRow=  TB[self.touchVedgeCol] + TB[self.touchLedgeCol] + TB[self.touchBedgeCol]



        return  TB[selectRow>0]





    def pureVclip(self,TB):

        """

        only keep those clouds that are clipped in the v axis, but not clipped in the l and b axis

        :param TB:
        :return:
        """

        selectTB= TB[ TB[self.touchVedgeCol] >0 ]

        selectTB= selectTB[ selectTB[self.touchLedgeCol] < 1  ]
        selectTB= selectTB[ selectTB[self.touchBedgeCol] < 1  ]

        return selectTB

    def pureLclip(self,TB):

        """

        only keep those clouds that are clipped in the v axis, but not clipped in the l and b axis

        :param TB:
        :return:
        """

        selectTB= TB[ TB[self.touchLedgeCol] >0 ]

        selectTB= selectTB[ selectTB[self.touchVedgeCol] < 1  ]
        selectTB= selectTB[ selectTB[self.touchBedgeCol] < 1  ]

        return selectTB

    def pureBclip(self,TB):

        """

        only keep those clouds that are clipped in the v axis, but not clipped in the l and b axis

        :param TB:
        :return:
        """

        selectTB= TB[ TB[self.touchBedgeCol] >0 ]

        selectTB= selectTB[ selectTB[self.touchLedgeCol] < 1  ]
        selectTB= selectTB[ selectTB[self.touchVedgeCol] < 1  ]

        return selectTB


    def pureLBclip(self,TB):

        """

        only keep those clouds that are not clipped in the v axis, but   clipped in the l or b axis

        :param TB:
        :return:
        """

        selectCriteria1 = TB[self.touchVedgeCol] < 1 #not clip in v axis

        selectCriteria2 = np.logical_or( TB[self.touchLedgeCol] > 0,    TB[self.touchBedgeCol] > 0,  )

        selectCriteria  = np.logical_and( selectCriteria1  ,    selectCriteria2  )

        selectTB= TB[  selectCriteria ]



        return selectTB



    def bothVAndLBClip(self,TB):
        """

        :return:
        """

        selectCriteria1 = TB[self.touchVedgeCol] > 0 #not in v axis

        selectCriteria2 = np.logical_or( TB[self.touchLedgeCol] >0,    TB[self.touchBedgeCol] > 0,  )

        selectCriteria  = np.logical_and( selectCriteria1  ,    selectCriteria2  )

        selectTB= TB[  selectCriteria ]



        return selectTB


    def noClipClouds(self,TB):

        selectTB= TB[ TB[self.touchVedgeCol]  < 1 ]

        selectTB= selectTB[ selectTB[self.touchLedgeCol] < 1  ]
        selectTB= selectTB[ selectTB[self.touchBedgeCol] < 1  ]

        return selectTB




    def drawFillingRelation(self,  fillingTB, calCode=None, drawCode="area"):
        """

        :param fillingTB:
        :return:
        """
        if calCode!=None:
            self.calCode=calCode

        ffTB= Table.read( fillingTB )
        ffTB=ffTB[ffTB[self.ffMWISPCol]>0]



        ffTB=self.addMWISPFFerror(ffTB)
        ffTB=self.addCfaFFerror(ffTB)


        #
        #ffTB=ffTB[ffTB[self.ffMWISPCol]>0]

        #ffTB=self.getCloudWithGoodFF(ffTB)

        #split into several samples


        #


        #only keep local clouds


        #ffTB=ffTB[abs(ffTB["y_cen"])>=2]


        print len(ffTB),"Number of total molecular clouds"

        pureVclipTB=self.pureVclip(ffTB)

        pureLBclipTB = self.pureLBclip( ffTB)

        noClipTB= self.noClipClouds(ffTB)
        bothClipTB = self.bothVAndLBClip(ffTB)

        edgeClodus= self.getEdgeClouds(ffTB)

        print "Number of edge clouds", len(edgeClodus)


        print  len(ffTB)- len(noClipTB ) -  len(pureLBclipTB )-  len(pureVclipTB ) - len(bothClipTB) ,"??????"


        fig = plt.figure(figsize=(10, 8))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})


        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        axFF = fig.add_subplot(1,1, 1)


        minVel = np.min(ffTB["v_cen"])
        maxVel = np.max(ffTB["v_cen"])

        #axFF.scatter(   ffTB["area_exact"]/3600., ffTB[self.ffMWISPCol]  , color="blue",  s=3)
        # ffMWISPCol   ffCfACol

        elinewidth = 0.6
        markerSize=2.1
        if drawCode==self.drawCodeArea: #angular aera

            self.drawErrorBar(axFF,noClipTB, drawCode =drawCode,markerSize=markerSize,color='gray',markerType=".",elinewidth=elinewidth,label="Complete in PPV space" ,showYError=False)
            self.drawErrorBar(axFF,pureVclipTB,drawCode =drawCode, markerSize=markerSize+0.8,color='b',markerType="D",elinewidth=elinewidth,label="Incomplete in v space" ,showYError=False)
            self.drawErrorBar(axFF,pureLBclipTB, drawCode =drawCode, markerSize=markerSize+0.8,color='r',markerType="^",elinewidth=elinewidth,label="Incomplete in l-b space" ,showYError=False)


            axFF.set_xlim([0, 100 ])
            axFF.set_xlabel("Angular area (square arcmin)")

        if drawCode==self.drawCodeSize: #angular size


            sc=self.drawErrorBar(axFF,noClipTB, drawCode =drawCode, markerSize=markerSize,color='gray',markerType=".",elinewidth=elinewidth,label="Complete in PPV space" ,showYError=False)
            #self.drawErrorBar(axFF,pureVclipTB, drawCode =drawCode, markerSize=markerSize+0.8,color='b',markerType="D",elinewidth=elinewidth,label="Incomplete in v space" ,showYError=False)
            #self.drawErrorBar(axFF,pureLBclipTB, drawCode =drawCode, markerSize=markerSize+0.8,color='r',markerType="^",elinewidth=elinewidth,label="Incomplete in l-b space" ,showYError=False)

            axFF.set_xlim([0,30 ])
            axFF.set_xlabel("Angular size (arcmin)")

            if 1:
                from mpl_toolkits.axes_grid1 import make_axes_locatable
                divider = make_axes_locatable(axFF)
                cax1 = divider.append_axes("right", size="3%", pad=0.05)
                cb = fig.colorbar(sc, cax=cax1)
                cax1.set_ylabel(r"$V_{\rm LSR}$ ($\rm km\ s^{-1}$)")





        if drawCode==self.drawCodeFlux: #angular size
            self.drawErrorBar(axFF,noClipTB, drawCode =drawCode, markerSize=markerSize,color='gray',markerType=".",elinewidth=elinewidth,label="Complete in PPV space" ,showYError=False)
            self.drawErrorBar(axFF,pureVclipTB, drawCode =drawCode, markerSize=markerSize+0.8,color='b',markerType="D",elinewidth=elinewidth,label="Incomplete in v space" ,showYError=False)
            self.drawErrorBar(axFF,pureLBclipTB, drawCode =drawCode, markerSize=markerSize+0.8,color='r',markerType="^",elinewidth=elinewidth,label="Incomplete in l-b space" ,showYError=False)



            axFF.set_xlim([0, 300 ])
            #axFF.set_xlabel("Angular size (arcmin)")
            axFF.set_xlabel(r"Flux ($\rm K\ km\ s$$^{-1}$ $\mathrm{\Omega}_\mathrm{A}$)")


        axFF.legend(loc=4)


        #### add Region to the figure


        at = AnchoredText(self.calCode , loc=1, frameon=False)
        axFF.add_artist(at)


        axFF.set_ylim([-0.2, 1.2  ])
        axFF.set_ylabel("Filling factor")

        plt.savefig( self.ffFigurePath+"{}FFindividual_{}.png".format(self.calCode, drawCode ), bbox_inches='tight', dpi=600)


    def checkRMS(self):
        """

        :return:
        """
        rms02FITS= "/home/qzyan/WORK/dataDisk/G2650/mosaic_U_rms_local.fits"
        rms016FITS= "/home/qzyan/WORK/dataDisk/G2650/mosaic_U_rms_localns.fits"

        #need to crop fits

        lRange=[34.5, 47.5]
        bRange = [-4.9, 4.9 ]

        rms02CropFITS= "rms02Crop.fits"
        rms016CropFITS= "rms016Crop.fits"

        doFITS.cropFITS2D(rms02FITS,outFITS=rms02CropFITS,Lrange=lRange,Brange=bRange, overWrite=True)
        doFITS.cropFITS2D(rms016FITS,outFITS=rms016CropFITS,Lrange=lRange,Brange=bRange, overWrite=True)



        rmsDataSm, rmsHeadSm=doFITS.readFITS( rms02CropFITS )

        rmsDataRaw, rmsheadRaw= doFITS.readFITS( rms016CropFITS )

        selectCriteria= np.logical_and(  rmsDataRaw > 0, rmsDataSm > 0 )


        goodSmValues = rmsDataSm[selectCriteria]
        goodRawValues= rmsDataRaw[selectCriteria]

        ratios= goodSmValues/goodRawValues

        print "Means ratio", np.mean( ratios )
        fits.writeto("rmsRatios.fits",rmsDataSm/rmsDataRaw,header=rmsHeadSm,overwrite=True )

        fig = plt.figure(figsize=(16, 9))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 16, 'serif': ['Helvetica']})
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        axScatter = fig.add_subplot(1,2, 1)
        axScatter.scatter(goodRawValues,goodSmValues,s=2,c='b')
        axScatter.set_xlabel( r"raw rms (0.16 $\rm km\ s^{-1}$)")
        #ax4.set_xlabel( r"$V_{\rm LSR}$ ($\rm km\ s^{-1}$)")

        axScatter.set_ylabel( r"smooth rms (0.2 $\rm km\ s^{-1}$)")

        #plot a

        theoryRatio=np.sqrt(0.16/0.2)

        xArray=np.asarray( [0.2,0.8] )



        axScatter.plot(xArray, xArray*theoryRatio,color='red')


        axHist = fig.add_subplot(1,2, 2)

        bins=np.arange(0.1,0.8,0.01)

        axHist.hist( goodRawValues,bins=bins, alpha=0.5,label="0.16, mean:{:.3f} K".format(np.mean(goodRawValues) ) )
        axHist.hist( goodSmValues,bins=bins, alpha=0.5,label="0.2, mean:{:.3f} K".format(np.mean(goodSmValues) ) )
        axHist.legend( loc=1 )

        axHist.set_xlim(0,1)

        print np.mean(goodSmValues) /np.mean(goodRawValues)


        plt.savefig("checkRMS.png" , bbox_inches='tight', dpi=600)


    def checkData(self):
        """
        draw several 12CO data, to compare the difference of spectral shape
        :return:
        """


        dataSm,headSm=doFITS.readFITS(self.localCO12FITS) #0.2 /km/s
        dataRaw,headRaw=doFITS.readFITS(self.rawLocalCO12FITS) #0.16 /km/s #replace this

        rmsSm, rmsHeadSm=doFITS.readFITS("/home/qzyan/WORK/dataDisk/G2650/mosaic_U_rms_local.fits")

        rmsRaw, rmsheadRaw= doFITS.readFITS("/home/qzyan/WORK/dataDisk/G2650/mosaic_U_rms_localns.fits")

        #draw four spectra lines

        fig = plt.figure(figsize=(10, 12))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 16, 'serif': ['Helvetica']})
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        ax1 = fig.add_subplot(4,1, 1)
        #axSm1 = fig.add_subplot(4,2, 2, sharex=axRaw1)

        #draw first spectra

        l,b=[37.9407 , 0.2716]
        self.drawSpectra(ax1,ax1,dataRaw,headRaw,dataSm,headSm,l,b)



        self.showNoise(ax1,rmsRaw,rmsheadRaw,rmsSm,rmsHeadSm,l,b)


        ax2 = fig.add_subplot(4,1, 2, sharex=ax1)
        #axSm2 = fig.add_subplot(4,2, 4, sharex=axRaw1)

        l,b=[27.9257170,-0.7087699]
        self.drawSpectra(ax2,ax2,dataRaw,headRaw,dataSm,headSm,l,b)
        self.showNoise(ax2,rmsRaw,rmsheadRaw,rmsSm,rmsHeadSm,l,b)

        ax3 = fig.add_subplot(4,1, 3, sharex=ax1)
        #ax3 = fig.add_subplot(4,1, 3, sharex=ax1 )

        l,b=[31.7994942,3.3178279]
        self.drawSpectra(ax3,ax3,dataRaw,headRaw,dataSm,headSm,l,b)

        ax3.legend(loc=5)
        self.showNoise(ax3,rmsRaw,rmsheadRaw,rmsSm,rmsHeadSm,l,b)





        ax4 = fig.add_subplot(4,1, 4, sharex=ax1)
        #ax4 = fig.add_subplot(4,1, 4, sharex=ax1)

        l,b=[48.1125907,0.3694292]
        self.drawSpectra(ax4,ax4,dataRaw,headRaw,dataSm,headSm,l,b)
        self.showNoise(ax4,rmsRaw,rmsheadRaw,rmsSm,rmsHeadSm,l,b)

        ax4.set_xlabel( r"$V_{\rm LSR}$ ($\rm km\ s^{-1}$)")
        ax4.set_xlabel( r"$V_{\rm LSR}$ ($\rm km\ s^{-1}$)")




        plt.savefig("checkData.png" , bbox_inches='tight', dpi=600)

    def showNoise(self,ax,rmsDataRaw,rmsHeadRaw,rmsDataSm,rmsheadSm,l,b):
        """

        :param ax:
        :param rmsDataRaw:
        :param rmsHeadRaw:
        :param rmsDataSm:
        :param rmsheadSm:
        :param l:
        :param b:
        :return:
        """

        noiseRaw= doFITS.getPixValue(rmsDataRaw,rmsHeadRaw,[l,b])
        noiseSm= doFITS.getPixValue(rmsDataSm,rmsheadSm,[l,b])

        ratio=  noiseSm/noiseRaw



        noiseStr="rms(0.16): {:.2f} K, rms(0.2): {:.2f} K, ratio:{:.3f}  (SQRT(0.16/0.2)=0.894)".format(noiseRaw,noiseSm,ratio)

        at = AnchoredText(noiseStr , loc=2, frameon=False)
        ax.add_artist(at)



    def drawSpectra(self,axRaw,axSm,dataRaw,headRaw,dataSm,headSm,l,b):


        sp1,v1=doFITS.getSpectraByLB( dataRaw, headRaw, l,b   )
        axSm.step(v1,sp1,color='green',lw=0.8 ,label=r"Raw (0.16 km s$^{-1}$)")

        #calculateNoise

        cutV=15
        spCut=sp1[v1>cutV]


        sp1,v1=doFITS.getSpectraByLB( dataSm, headSm, l,b   )
        axRaw.step(v1,sp1,color='blue',lw=0.8, label=r"Smooth (0.2 km s$^{-1}$)")
        spCut=sp1[v1>cutV]

        #
        lbStr=r"(l,b)=({:.2f}$^\circ$,{:.2f}$^\circ$)".format(l,b)
        at = AnchoredText( lbStr , loc=4, frameon=False)
        axRaw.add_artist(at)
        axRaw.set_ylabel("T$_{\rm mb}$ (K)")
        #at = AnchoredText( lbStr , loc=1, frameon=False)
        #axSm.add_artist(at)
        #axSm.set_ylabel("T$_{\rm mb}$ (K)")




    def selectBySizeAndFactorRange(self,TBFileName,sizeRange=None,factorRange=None):
        """
        Used to select sub samples
        :param sizeRange:
        :param factorRange:
        :return:
        """
        TB=Table.read( TBFileName )

        sizes= self.getCloudSize( TB )

        if sizeRange!=None:
            TB=self.selectTBByColRange(TB,sizes,sizeRange)

        if factorRange !=None:


            TB=self.selectTBByColRange(TB,TB[self.ffMWISPCol],factorRange)

        return TB


    def selectTBByColRange(self,TB,col,colRange):
        """
        used to select table by columan range
        :param TB:
        :param col:
        :param colRange:
        :return:
        """

        selectCriteria1 = col >= np.min(colRange)
        selectCriteria2 = col <= np.max(colRange)

        selectCriteria= np.logical_and(selectCriteria1,selectCriteria2)

        return TB[selectCriteria]




    def getCloudSize(self,TB):

        return np.sqrt(   TB["area_exact"]/np.pi )*2

    def drawErrorBar(self,ax, TB,drawCode="area", showYError=False, color="gray", markerType='.', markerSize=2,elinewidth=0.6,label=""):


        if drawCode== self.drawCodeArea:
            

            drawX= TB["area_exact"]

        if drawCode == self.drawCodeSize:


            drawX=np.sqrt(   TB["area_exact"]/np.pi )*2

        if drawCode == self.drawCodeFlux:
            drawX= TB["sum"]*self.getVelResolution()


        #
        drawY=  TB[self.ffMWISPCol]
        #drawY1,error1=self.getFillingErrorAndError(5,TB)
        #drawY2,error2=self.getFillingErrorAndError(1,TB)
        #drawY= drawY1/drawY2

        yerr=  TB[self.ffMWISPErrorCol]

        if not showYError:
            yerr =  None
        import matplotlib.cm as cm
        import matplotlib as mpl

        Vs=  TB["v_cen"]

        cmap = plt.cm.jet
        # norm = matplotlib.colors.BoundaryNorm(np.arange(0,30,0.1), cmap.N)
        normV = mpl.colors.Normalize(vmin=np.min(Vs), vmax=np.max(Vs) )
        #normV = mpl.colors.Normalize(vmin=np.min(Vs), vmax=np.max(Vs) )


        m = plt.cm.ScalarMappable(norm=normV, cmap=cmap)
        v_color = np.array([(m.to_rgba(v)) for v in Vs])

        if color=="gray":
            #c =Vs, cmap=cmap,norm=normV,
            #ax.errorbar(drawX,   drawY , yerr= yerr ,    markersize=markerSize , linestyle='none', c= v_color ,     marker= markerType , capsize=0,  elinewidth=elinewidth , lw=1, label= label )
            #ax.errorbar(drawX,   drawY , yerr= yerr ,    markersize=markerSize , linestyle='none', c= Vs ,  cmap=cmap,norm=normV,   marker= markerType , capsize=0,  elinewidth=elinewidth , lw=1, label= label )
            sc = ax.scatter(drawX, drawY, c=Vs, cmap=cmap, norm=normV, s=3, facecolors='none', lw=0.5, marker="o",label="")
            #sc = ax.scatter(drawX, drawY, c=color, cmap=cmap, norm=normV, s=3, facecolors='none', lw=0.5, marker="o",label="")

            return sc



        else:

            ax.errorbar(drawX,   drawY , yerr= yerr ,    markersize=markerSize , linestyle='none', c= color ,    marker= markerType , capsize=0,  elinewidth=elinewidth , lw=1, label= label )


    def produceLVByRemovingTheLargestCloud(self):

        """
        remove the largest molecular cloud, to see, if there is any regular variation in the lv Diagram
        :return:
        """

        rawFITS=self.getRawCOFITS( )

        cleanFITS = self.getSmoothAndNoiseFITSSingle(getCleanFITS=True)

        removeIDs=[146]


        labelData,labelHead=doFITS.readFITS(cleanFITS)

        rawData,rewhead= doFITS.readFITS(rawFITS)



        for eachID in removeIDs:
            labelData[labelData==eachID]=0 #remove those clouds

        rawData[labelData==0]=0

        pvData,pvHead=doFITS.readFITS("/home/qzyan/WORK/diskMWISP/fillingFactorData/data/rawLocalCO12_LV.fits")

        lvData=np.sum(rawData,axis=1)



        print rawFITS
        print cleanFITS


        fits.writeto("SmallCloudPV.fits",lvData, header=pvHead,overwrite=True)



    def drawIntOveralIntMap(self,fitsFile, showLineStr ):
        fig = plt.figure(figsize=(10, 6))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 13, 'serif': ['Helvetica']})
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]



        dataRaw, headRaw= doFITS.readFITS(fitsFile)
        
        wcsCO=WCS(headRaw)
        axCO = pywcsgrid2.subplot(111, header=wcsCO)


        
        #dataLabel,headLabel=
        Nz,Ny,Nx= dataRaw.shape

        noiseSingleChannel=self.getMeanRMS()
        velRes=self.getVelResolution()

        intData=np.nansum(dataRaw,axis=0)*velRes

        guessNoise=np.sqrt(Nz)*velRes*noiseSingleChannel


        cmapCO = plt.cm.bone
        cmapCO.set_bad('black')
        
        axCO.imshow(  intData , origin='lower', cmap=cmapCO, vmin=guessNoise, vmax=guessNoise*6, interpolation='none')
        axCO.set_ticklabel_type("absdeg", "absdeg")
        axCO.axis[:].major_ticks.set_color("w")

        #get spectra

        spec0,vel0=doFITS.getSpectraByIndex(dataRaw,headRaw,0,0)
        v0=vel0[0]
        v1=vel0[-1]

        vInfo= "Integrate from {:.2f} to {:.2f} km/s".format(v0, v1)
        at = AnchoredText( vInfo, loc=1, frameon=False ,prop={"color":"white"} )
        axCO.add_artist(at)
        at = AnchoredText( showLineStr , loc=2, frameon=False, prop={"color": "white"})
        axCO.add_artist(at)
        plt.savefig(self.calCode+"OverAllIntMap.png" , bbox_inches='tight', dpi=600)


    def drawSingleChannelMap(self,fitsFile,channelN,showLineStr):
        plt.clf()

        fig = plt.figure(figsize=(10, 6))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 13, 'serif': ['Helvetica']})
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]




        dataRaw, headRaw = doFITS.readFITS(fitsFile)

        wcsCO = WCS(headRaw)
        axCO = pywcsgrid2.subplot(111, header=wcsCO)

        # dataLabel,headLabel=
        Nz, Ny, Nx = dataRaw.shape

        noiseSingleChannel = self.getMeanRMS()
        velRes = self.getVelResolution()

        intData = dataRaw[channelN]* velRes  #np.nansum(dataRaw, axis=0) * velRes

        guessNoise =  velRes * noiseSingleChannel

        cmapCO = plt.cm.bone
        cmapCO.set_bad('black')

        axCO.imshow(intData, origin='lower', cmap=cmapCO, vmin=guessNoise, vmax=guessNoise * 6, interpolation='none')
        axCO.set_ticklabel_type("absdeg", "absdeg")
        axCO.axis[:].major_ticks.set_color("w")

        # get spectra

        spec0, vel0 = doFITS.getSpectraByIndex(dataRaw, headRaw, 0, 0)
        v0 = vel0[channelN]


        vInfo = "VLSR = {:.2f} km/s".format(v0 )
        at = AnchoredText(vInfo, loc=1, frameon=False, prop={"color": "white"})
        axCO.add_artist(at)

        at = AnchoredText( showLineStr , loc=2, frameon=False, prop={"color": "white"})
        axCO.add_artist(at)


        plt.savefig(self.calCode + "channelMap{}.png".format(channelN),  bbox_inches='tight', dpi=600)

    def testBadChannel(self,showChannel ,lRange,bRange):
        #self.calCode=self.codeRawScuCO13

        #lRange,bRange=doFITS.box(41.9992119,2.9129983,1724.880 ,1540.763 ,0)

        rawFITS=self.getRawCOFITS()

        dataRaw,headRaw= doFITS.readFITS(rawFITS)

        spec,vel=doFITS.getAverageSpecByLBrange(rawFITS, lRange, bRange )
        
        #spec,vel=doFITS.getSpectraByLB(dataRaw,headRaw, 42.0832649, 2.8251944 )

        lineC18OStr= r"$\mathrm{C}^{18}\mathrm{O}(J=1\rightarrow0)$"


        showLine=lineC18OStr

        plt.clf()

        fig = plt.figure(figsize=(10, 8))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 13, 'serif': ['Helvetica']})
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]


        wcsCO=WCS(headRaw)

        #showChannel=165

        axCO = pywcsgrid2.subplot(211, header=wcsCO)

        cmapCO = plt.cm.bone
        cmapCO.set_bad('black')
        axCO.imshow( dataRaw[showChannel], origin='lower', cmap=cmapCO, vmin=0, vmax=3, interpolation='none')

        axCO["gal"].plot([lRange[0], lRange[0]], bRange, 'r--', lw=0.5)
        axCO["gal"].plot([lRange[1], lRange[1]], bRange, 'r--', lw=0.5)

        axCO["gal"].plot(lRange, [bRange[0], bRange[0]], 'r--', lw=0.5)
        axCO["gal"].plot(lRange, [bRange[1], bRange[1]], 'r--', lw=0.5)

        ###################
        ExtendDeg = 2 # 1.5

        cutlSmall = min(lRange) - ExtendDeg
        cutlLarge = max(lRange) + ExtendDeg

        cutbSmall = min(bRange) - ExtendDeg
        cutbLarge = max(bRange) + ExtendDeg

        backWCS = WCS(headRaw, naxis=2)  # WCS(allHead)
        x0, y0 = backWCS.all_world2pix(cutlLarge, cutbSmall, 0)
        #
        x1, y1 = backWCS.all_world2pix(cutlSmall, cutbLarge, 0)

        xRange = [x0, x1];
        yRange = [y0, y1]
        axCO.set_xlim(min(xRange), max(xRange)),
        axCO.set_ylim(min(yRange), max(yRange)),
        #################
        axSpectra = fig.add_subplot(2,1, 2)

        axCO.set_ticklabel_type("absdeg", "absdeg")
        axCO.axis[:].major_ticks.set_color("w")

        axSpectra.step( range(len(spec)),spec,lw=0.5,color='b' ,where="mid")

        axSpectra.set_xlabel("Channel number")
        axSpectra.set_ylabel("Brightness temperature (K)")

        spec0, vel0 = doFITS.getSpectraByIndex(dataRaw, headRaw, 0, 0)
        v0 = vel0[showChannel]
        vInfo = "{:.2f} km/s ".format(v0,showChannel )
        at = AnchoredText(vInfo, loc=1, frameon=False, prop={"color": "white"})
        axCO.add_artist(at)

        vInfo = "Channel Number {}".format( showChannel )
        at = AnchoredText(vInfo, loc=4, frameon=False, prop={"color": "white"})
        axCO.add_artist(at)

        axSpectra.axvline(x=showChannel, ls="--", color='black',lw=0.5,alpha=0.5)

        
        at = AnchoredText( showLine , loc=1, frameon=False, prop={"color": "black"})
        axSpectra.add_artist(at)


        at = AnchoredText( showLine , loc=2, frameon=False, prop={"color": "white"})
        axCO.add_artist(at)


        plt.savefig("testBadChannel.png" , bbox_inches='tight', dpi=600)


        plt.clf()
        self.drawIntOveralIntMap(rawFITS,showLine)
        self.drawSingleChannelMap(rawFITS,showChannel,showLine)




    def getConsecutiveMask(self,dataCO,rmsData,cutoff=3,conNumber=3):
        """
        return a mask that above cutoff*3, sigma

        find the points that contains three consecutivechannels largeer than 3 sigma
        :param dataCO:
        :param rmsData:
        :return:
        """

        Nz,Ny,Nx=dataCO.shape

        zero2D=np.zeros((1,Ny,Nx))
        sigmaData=dataCO/rmsData

        cutoffmask= sigmaData  >=cutoff

        extendedArray=np.vstack( (zero2D,cutoffmask,zero2D) )

        sumArray = extendedArray[0:-2,:,:] +  extendedArray[1:-1,:,:]
        sumArray = sumArray + extendedArray[2: ,:,:]


        consecutiveData = sumArray>=conNumber #consecutively 3 channnel l



        return consecutiveData*1


    def removeBadSgrChannelCO12(self):
        """

        :return:
        """
        ##
        rawCOFITS="/home/qzyan/WORK/diskMWISP/fillingFactorData/data/rawSgrCO12BackUp.fits"
        labelFITS="/home/qzyan/WORK/diskMWISP/fillingFactorData/data/rawSgrCO12rawSgrCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_CleanBackUp.fits"

        saveCorrectedFITS = "/home/qzyan/WORK/diskMWISP/fillingFactorData/data/rawSgrCO12Correct.fits"

        rawCOData,rawCOHead= doFITS.readFITS( rawCOFITS)
        labelData, rawCOHead= doFITS.readFITS( labelFITS)

        self.calCode=self.codeRawSgrCO12

        rmsFITS=self.getRMSFITS()
        #27,28,29,three bad channels

        dataNoise=self.produceNoiseFITS(3,rmsFITS=rmsFITS)
        
        print dataNoise.shape

        channel26Label = labelData[26]
        channel30Label = labelData[30]

        #positions, need to mask
        maskLabel =  np.logical_and(channel26Label==0, channel30Label==0  )

        rawCOData[27][maskLabel]=dataNoise[0][maskLabel]
        rawCOData[28][maskLabel]=dataNoise[1][maskLabel]
        rawCOData[29][maskLabel]=dataNoise[2][maskLabel]


        fits.writeto(saveCorrectedFITS,rawCOData,header=rawCOHead,overwrite=True )



    def calculateChiSquare(self,x,y,yError, ffFunction, parameters):
        """

        :param x:
        :param y:
        :param yError:
        :param function:
        :param parameters:
        :return:
        """


        yMu= ffFunction(x, *parameters )

        x2=(yMu-y)**2/yError**2

        x2=np.sum(x2)

        print x2,"chisquare"

        return x2


    def weightedRMS(self,x,y,yError, ffFunction, parameters):
        """

        :param x:
        :param y:
        :param yError:
        :param function:
        :param parameters:
        :return:
        """


        yMu= ffFunction(x, *parameters )


        residual= yMu-y
        #x2=(yMu-y)**2/yError**2

        #x2=np.sum(x2)


        weights= 1./yError**2
        weights=weights/np.sum( weights )



        return np.sqrt( np.sum(weights*residual**2) )



    def testThreeFunctions(self,testFun,fillingTB,errorCut=None):
        ##calculate the chi-square

        ffTB=self.addFFValues(fillingTB)

        if errorCut is not None:

            selectByError=  ffTB[self.ffMWISPErrorCol] /ffTB[self.ffMWISPCol]
            ffTB=ffTB[selectByError<=errorCut]
        else:
            pass
        drawX = self.getCloudSize(ffTB)
        drawY = ffTB[self.ffMWISPCol]
        yError= ffTB[self.ffMWISPErrorCol]


        para,paraError=self.fittingFFAndSize(testFun, drawX,drawY, yError,sizeError=None,  useODR=False)
        print para,paraError
        return self.calculateChiSquare(drawX,drawY,yError, testFun,para)


    def splitTBIntoTrainingTest(self,TB,trainingRatio=0.8):
        """
        :param TB:
        :return:
        """
        N=len(TB)

        #
        indexArray = np.arange(N)


        np.random.shuffle(indexArray)

        cutIndex=int(round( N/1.0*trainingRatio ) )


        newRandonTB= TB[indexArray]

        return newRandonTB[:cutIndex] , newRandonTB[cutIndex:]




    def testThreeFunctionsByTraining(self, testFun, fillingTB,errorCut=None , trainingRatio=0.7,drwaFigure=False):
        ##calculate the chi-square

        ffTB = self.addFFValues(fillingTB)


        if errorCut is not None:

            selectByError = ffTB[self.ffMWISPErrorCol] / ffTB[self.ffMWISPCol]
            ffTB = ffTB[selectByError <= errorCut]
        else:
            pass

        #split the table with training and test

        trainingTB,testTB= self.splitTBIntoTrainingTest(ffTB,trainingRatio=trainingRatio)



        drawX = self.getCloudSize(trainingTB)
        drawY = trainingTB[self.ffMWISPCol]
        yError = trainingTB[self.ffMWISPErrorCol]

        para, paraError = self.fittingFFAndSize(testFun, drawX, drawY, yError, sizeError=None, useODR=False)
        print para, paraError

        testX= self.getCloudSize(testTB)
        testY=  testTB[self.ffMWISPCol]
        testYerror=  testTB[self.ffMWISPErrorCol]

        saveName="paraN{}ErrorCut{}TrainingRatio{}".format(len(para),errorCut,trainingRatio)
        
        if drwaFigure: #draw a figure
            plt.clf()
            fig = plt.figure(figsize=(10, 8))
            rc('text', usetex=True)
            rc('font', **{'family': 'sans-serif', 'size': 13, 'serif': ['Helvetica']})
            mpl.rcParams['text.latex.preamble'] = [
                r'\usepackage{tgheros}',  # helvetica font
                r'\usepackage{sansmath}',  # math-font matching  helvetica
                r'\sansmath'  # actually tell tex to use it!
                r'\usepackage{siunitx}',  # micro symbols
                r'\sisetup{detect-all}',  # force siunitx to use the fonts
            ]

            axCO = plt.subplot(111)

            axCO.scatter(drawX,drawY,color="gray",s=5  )

            axCO.scatter(testX,testY,color="blue",s=5  )
            axCO.set_xlim(0,30)
            axCO.legend(loc=4)

            drawRange=  np.linspace(0,30,100)

            axCO.plot( drawRange, testFun(drawRange,*para),"r-",lw=1 )


            axCO.set_xlabel("Error threshold")
            axCO.set_ylabel("Chi-Square")

            plt.savefig(saveName+"trainingTest.png", bbox_inches='tight', dpi=300)

        #return self.calculateChiSquare(testX, testY, testYerror, testFun, para)
        return self.weightedRMS(testX, testY, testYerror, testFun, para)

    def getLabelSet(self,labelData):
        """

        :param labelData:
        :return:
        """
        noiseMask= np.min(   labelData[0] )

        clusterIndex1D = np.where(labelData > noiseMask )
        clusterValue1D = labelData[clusterIndex1D]

        Z0, Y0, X0 = clusterIndex1D
        labelSets=[Z0, Y0, X0, clusterValue1D ]

        return  labelSets

    def getSmoothFluxColCutoff(self, rawCOFITS, labelFITS, rmsFITS,TBFile ):
        """
        get the change of flux according to the cutoff sigmas, remember, this has to be cut according to the rmsData
        :return:
        """

        ####
        print "Extracting cutoff flux from ", rawCOFITS

        dataLabel , headLabel = doFITS.readFITS( labelFITS )
        dataRaw , headRaw = doFITS.readFITS( rawCOFITS )
        rmsData,rmsHead=doFITS.readFITS( rmsFITS )

        TB=Table.read(TBFile)
        self.addCutOffFluxColnames(TB)


        labelSets= self.getLabelSet(dataLabel) #[Z0, Y0, X0, clusterValue1D ]



        #colName,colPixName=self.getFluxColNameNoiseChange(cutoffFactor)

        widgets = ['Extracting cutoff flux:', Percentage(), ' ', Bar(marker='0', left='[', right=']'), ' ', ETA(), ' ',
                   FileTransferSpeed()]  # see docs for other options
        pbar = ProgressBar(widgets=widgets, maxval=len(TB))
        pbar.start()

        TB.sort("area_exact")
        TB.reverse()
        i=0
        for eachRow in TB:
            i=i+1
            #gc.collect()

            ID = eachRow["_idx"]

            cloudIndex,cloudIndex2D = self.getIndices(labelSets, ID,return2D=True)  # np.where( cleanDataSM1==ID )
            coValues=   dataRaw[cloudIndex]
            rmsValues= rmsData[cloudIndex2D]


            for eachCutoff in self.cutoffFactors:

                selectCoValues= coValues[ coValues>=rmsValues*eachCutoff ]

                #print selectCoValues

                colName, colPixName = self.getFluxColNameCutoff(eachCutoff)

                eachRow[colName]= np.sum( selectCoValues )*self.getVelResolution()
                eachRow[colPixName]= len(selectCoValues)


            



            pbar.update(i)

        pbar.finish()
        TB.write("cutoffFlux_"+TBFile , overwrite=True)

        TB=self.calculateFillingFactorCutoff(TB)
        TB.write("cutoffFF_"+TBFile , overwrite=True)

        print TB
        #now TB contains all  cutOffFlux Info to calculate fifling factors


    def cleanRawData(self):

        """
        remove data with largeRMS
        :return:
        """

        #lRange, bRange, that need to set as nan

        #lRangeNAN=[2814, 2868] # 26.25
        #bRangeNAN=[1062,1122] #3.75 4.25
        for eachCode in self.allRawCodeList:

            self.calCode=eachCode


            rmsFITS= self.getRMSFITS()
            rmsData,rmsHead=doFITS.readFITS( rmsFITS )

            wcsRMS=WCS(rmsHead,naxis=2)

            lowerLeftPoint = wcsRMS.wcs_world2pix(26.25,3.75,0)
            upperRightPoint= wcsRMS.wcs_world2pix(25.8,4.25,0)


            lowerLeftPoint=map(round,lowerLeftPoint )
            upperRightPoint=map(round,upperRightPoint )


            x0,y0=map(int,lowerLeftPoint )
            x1,y1=map(int,upperRightPoint )

            rmsData[y0+1:y1, x0+1:x1+1 ]=np.nan
            fits.writeto(rmsFITS,rmsData,header=rmsHead,overwrite=True)

            COFITS=self.getRawCOFITS()

            dataCO,headCO=doFITS.readFITS(COFITS)
            dataCO[:,y0+1:y1, x0+1:x1+1 ]=np.nan

            fits.writeto( COFITS ,dataCO,header=headCO,overwrite=True)




    def zzz(self):
        """

        :return:
        """

