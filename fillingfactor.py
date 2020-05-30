
import numpy as np
import matplotlib as mpl
mpl.use('agg')

import radio_beam
from spectral_cube import SpectralCube
from astropy import units as u

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

sys.path.insert(1, '/home/qzyan/WORK/myDownloads/MWISPcloud')
sys.path.insert(1, '/home/yqzpmo/github/MWISPdbscan')

from onlyDBSCAN import allDBSCAN

doAllDBSCAN = allDBSCAN()

doFITS=myFITS() #used to deal fits with myPYTHON



########
def ffFunction(x,a, b, c):
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
def aaaa():
    pass



class checkFillingFactor(object):

    rootPath="./"

    #saveFITSPath="/home/qzyan/WORK/diskMWISP/fillingFactorData/"
    saveFITSPath="./"

    dataPath= saveFITSPath+"data/"
    rawFITS=  dataPath + "G2650Local30.fits"
    rmsFITSPath= rootPath+"rmsFITSpath/"
    #tmpPath= rootPath +"tmpFiles/"
    tmpPath= saveFITSPath +  "tmpFiles/"
    #########out
    figurePath=  rootPath +"figurePath/"

    intFigurePath = rootPath+"cloudIntMap/"

    cloudCubePath= saveFITSPath+ "cloudCubes/"
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

    allCodeList =  localCodeList+  sgrCodeList +  scuCodeList +  outCodeList
    allNoiseList = noiseList + noiseList + noiseList + noiseList
    allArmList=  ["local"]*3+  ["Sagittarius"]*3+  ["Scutum-Centaurus"]*3+  ["Outer"]*3    #armStrList+armStrList+armStrList+armStrList

    allLineList= lineStrList+ lineStrList+lineStrList+lineStrList

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


    #################################
    touchLedgeCol= "touchLedge"
    touchBedgeCol= "touchBedge"
    touchVedgeCol= "touchVedge"

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


        observedData=dataSM+noiseData


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







    def drawCloudNumberChange(self):

        testCode=self.codeLocalCO12
        self.calCode = testCode
        TB,TBFiles=self.getAbsNoiseCleanTBList( absK=self.MWISPrmsCO12)

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

        beamSize=self.smoothFactors*self.getBeamSize(self.calCode)

        #self.drawNoiseEffectSingle(axSens,1.5)

        #axBeam.plot(beamSize,Nlist,'o-',color='b')


        mwispF,cfaF,parameters=self.getFillingFactorAndDraw(  beamSize, Nlist, 1, drawFigure=False)
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





    def getBeamSize(self,calCode="3333"):
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

    def getCloudCubes(self, calCode,rawCOFITS, labelsFITS,cloudTBFile ,writeFITS=False ):

        """

        #output all data cubes for each cloud

        :return:
        """

        #################
        self.calCode=calCode
        savePath =  self.checkCloudCubeSavePath() #"./cloudSubCubes/"



        cloudTB = Table.read(cloudTBFile)
        dataCluster, headCluster = myFITS.readFITS(labelsFITS)
        dataCO, headCO = myFITS.readFITS( rawCOFITS)


        minV = np.nanmin(dataCluster[0])
        wcsCloud = WCS(headCluster)
        clusterIndex1D = np.where(dataCluster > minV)
        clusterValue1D = dataCluster[clusterIndex1D]
        Z0, Y0, X0 = clusterIndex1D

        fitsZero = np.zeros_like(dataCluster,dtype=np.float)


        #add a function that check if the cloud touches L edge,B edges, or zEdges

        Nz, Ny, Nx = dataCO.shape

        ######
        cloudTB[self.touchLedgeCol]=cloudTB["peak"]*0
        cloudTB[self.touchBedgeCol]=cloudTB["peak"]*0
        cloudTB[self.touchVedgeCol]=cloudTB["peak"]*0

        widgets = ['Extracting cloud cubes: ', Percentage(), ' ', Bar(marker='0', left='[', right=']'), ' ', ETA(), ' ',
                   FileTransferSpeed()]  # see docs for other options
        pbar = ProgressBar(widgets=widgets, maxval=len(cloudTB))
        pbar.start()
        i=0
        for eachC in cloudTB:

            i=i+1

            cloudID = eachC["_idx"]
            saveName = "{}cloud{}cube.fits".format(calCode,cloudID)

            cloudIndex = self.getIndicesRaw(Z0, Y0, X0, clusterValue1D, cloudID)
            fitsZero[cloudIndex] = dataCO[cloudIndex]

            cloudZ0, cloudY0, cloudX0 = cloudIndex

            minZ = np.min(cloudZ0)
            maxZ = np.max(cloudZ0)

            minY = np.min(cloudY0)
            maxY = np.max(cloudY0)

            minX = np.min(cloudX0)
            maxX = np.max(cloudX0)

            ####
            if minZ==0 or maxZ==Nz:
                eachC[self.touchVedgeCol]=1


            if minY==0 or maxY==Ny:
                eachC[self.touchBedgeCol]=1

            if minX==0 or maxX==Nx:
                eachC[self.touchLedgeCol]=1



            cropWCS = wcsCloud[minZ:maxZ + 1, minY:maxY + 1, minX:maxX + 1]

            cropData = fitsZero[minZ:maxZ + 1, minY:maxY + 1, minX:maxX + 1]

            saveFullName=os.path.join(savePath,saveName)

            if writeFITS:
                fits.writeto(saveFullName , cropData, header=cropWCS.to_header(), overwrite=True)

            fitsZero[cloudIndex] = 0

            pbar.update(i)
        pbar.finish()
        cloudTB.write("edgeInfo_"+cloudTBFile,overwrite=True)

    def getTotalFlux(self,fitsFile):
        """
        Use DBSCAN to maked to
        :param fitsFile:
        :return:
        """




    def smoothFITSbySMFactor(self,rawCOFITS, rawBeamSize = 52./60  ):

        """
        Default
        :param rawCOFITS:
        :param rawBeamSize:
        :return:
        """
        #calCode=rawCOFITS[0:-5]
        #self.calCode=calCode
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
            totalFlux=np.sum(eachTB["sum"])*0.2 #convert to 0.2 km/s

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

    def cleanFITSsigma2(self,FITSfile,cutoff=2,minPts=4,contype=1,removeFITS=False):
        """
        Used  DBSCAN to mask noise pixels, then calculate the total flux
        :param FITSfile:
        :return:
        """
        print "Cleaning FITS ",FITSfile
        saveTag=os.path.basename(FITSfile)
        saveTag = saveTag[0 : -5]

        #rawFITS="/home/qzyan/WORK/dataDisk/MWISP/G120/100_150_U.fits" #input raw coFITS
        #doAllDBSCAN.rmsCO12=0.5 # K, set the rms



        FITSrms =  self.getMeanRMS() #doFITS.getRMSFITS(FITSfile , ""  , returnRMSValue=True )

        doAllDBSCAN.rmsCO12=  FITSrms
        print "The rms ", FITSrms
        #saveTag="Q2Test" #your project code

        #cutoff=2 #in units of sigma
        #minPts=4
        #contype=1

        doAllDBSCAN.testAllPath=self.tmpPath
        #doAllDBSCAN.testAllPath=  FITSrms
        doAllDBSCAN.setfastDBSCANrms(FITSrms)

        dbscanLabelFITS, rawDBSCANTBFile  = doAllDBSCAN.pureDBSCAN(FITSfile, cutoff , MinPts=minPts, saveTag= saveTag , connectivity= contype , inputRMS=FITSrms, redo=True, keepFITSFile=True)


        cleanFITSlabel = doAllDBSCAN.getMaskByLabel(  dbscanLabelFITS  , rawDBSCANTBFile, onlyClean=True )

        TB=doAllDBSCAN.getMaskByLabel( dbscanLabelFITS ,  rawDBSCANTBFile , onlySelect=True ) #only Select=True, meas we only want to select the catlaog
        TB.write(  cleanFITSlabel[0:-5]+".fit"  , overwrite=True)
        if removeFITS:
            os.remove(cleanFITSlabel)
            os.remove(dbscanLabelFITS)



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



    def getSmoothAndNoiseFITSSingle( self, smFactor=1.0, noiseFactor=0.0, calCode=None, getCleanTBFile=False,dbscanCode="dbscanS2P4Con1" ,getCleanFITS=False):
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

                sumV= np.nansum(data[data >= downToK])*0.2

                sumList.append( sumV )

            return np.asarray( sumList )


        else:
            data,head=myFITS.readFITS(fitsName)
            self.maskNoiseEdge(data)
            #print np.nansum(data)
            return np.nansum(data[data>=downToK])*0.2


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

            dirrectCutValue = np.nansum(data[data>eachK])*0.2

            #if we take the dirrect cutValue as the 2 sigma
            noiseData=self.addNoiseSingle( smoothedFITS ,   absoluteNoise= eachK/4. , returnData=True)

            noiseDataSelect = noiseData[noiseData>eachK]

            print "Total Number of summed voxels, ", len( noiseDataSelect )

            noiseCutValue = np.nansum(noiseDataSelect)*0.2

            diff=   noiseCutValue - dirrectCutValue
            print  dirrectCutValue,  noiseCutValue,   diff/dirrectCutValue


    def fluxAboveCutoff(self,data,cutoff):
        """
        sum all flux above the cutoff threshold
        :param data:
        :param cutoff:
        :return:
        """

        return np.nansum(data[data>cutoff])*0.2



    def getEqualCutOff(self,fitsFile,TBFile):

        """
        get the equivalent cutoff of fitsFile to get the equivalent flux with DBSCAN catalg
        :param fitsFile:
        :param TBFile:
        :return:
        """
        print fitsFile
        print TBFile

        tbFile=Table.read(TBFile)

        dbscanFlux = np.sum(tbFile["sum"])*0.2

        precision=0.0001 #percent



        downToKUp = 3. #upper, the lower is 1 K
        downToKLow = 1. #upper, the lower is 1 K


        data, head = myFITS.readFITS(fitsFile)
        self.maskNoiseEdge(data)



        while 1:
            #calcualate cutoff rms

            downToKTry = (downToKUp + downToKLow) / 2.
            #print "Trying downToK, ", downToKTry
            sumV = np.nansum(data[data >= downToKTry]) * 0.2

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




    def getFillingFactorAndDraw(self,beamList,fluxList,calID,drawFigure=False):
        """

        :param beamList:
        :param fluxList:
        :return:
        """

        x=np.asarray( beamList )

        y = np.asarray( fluxList )


        #print x
        #print y
        try:
            params, paramas_covariance = optimize.curve_fit(ffFunction, x, y, p0=[np.mean(y), 0.5, np.mean(y)])
        except:

            return 0, 0, [[0,0,0],[0,0,0]]



        errors = np.sqrt(np.diag(paramas_covariance))

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
        axFitting.scatter(  x  ,  y , s=15,  color='red'   )

        fittingX = np.arange(0, np.max(x), 0.01)
        # axFitting.plot( fittingX  ,  ffFunction(fittingX,params[0], params[1] , params[2]   ), color='blue'  )
        axFitting.plot(fittingX, ffFunction(fittingX, *params), color='blue', lw=1.5)

        fillingWMISP = "The filling factor ({}) of MWISP is: {:.3f}".format(self.calCode, wmsipFilling)
        fillingCfA = "The filling factor ({}) of CfA is: {:.3f}".format(self.calCode, cfaFilling)

        # at = AnchoredText(r"The noise rms is {} K (1$\sigma$)".format(self.MWISPrmsCO12), loc=1, frameon=False)
        # axFitting.add_artist(at)

        # if cfaBeam> np.max(x):
        # cfaFilling= 0

        at = AnchoredText(fillingWMISP + "\n" + fillingCfA, loc=1, frameon=False)
        axFitting.add_artist(at)
        ###

        axFitting.axvline(x=0, ls="--", color='black')

        # draw CfA resolution
        axFitting.axvline(x=cfaBeam, ymax=0.2, ls="--", color='green', label="CfA resolutuion (8.5 arcmin)")
        axFitting.legend(loc=5)
        ###########################
        saveTag= "{}_factorFitting_ID{}".format(self.calCode,calID)

        plt.savefig(self.figurePath+"{}.png".format( saveTag ), bbox_inches='tight', dpi=200)

        plt.close(fig )
        gc.collect()
        return wmsipFilling, cfaFilling,  fittingParaAndError



    def drawFillingFactor(self, absK=0.5, useArea=False ):

        """
        Get clean TB files to draw total flux change
        :param calCode:
        :param absK:
        :return:
        """


        ##get clean file

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

        plt.savefig("{}.pdf".format( saveTag  ), bbox_inches='tight')
        plt.savefig("{}.png".format( saveTag ), bbox_inches='tight', dpi=600)
        #the flux is saved in npy
        return   wmsipFilling, cfaFilling  , fittingParaAndError


    def printFillingCat(self):
        """
        print the filling factor for each arm, each lines and together with the fitting parameters
        :return:
        """

        ###

        mwispFillingList=np.load("mwispFilling.npy")
        cfaFillingList=np.load("cfaFilling.npy")
        paraFitting = np.load("fittingPara.npy")


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
        mwispFillingList=np.load("mwispFilling.npy")
        cfaFillingList=np.load("cfaFilling.npy")
        paraFitting = np.load("fittingPara.npy")





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

            codeI=self.allCodeList[i]
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
            fluxList_1.append( np.sum(data)*0.2*omega  )

            data[data<2*rms]=0
            fluxList_2.append( np.sum(data)*0.2*omega  )

            data[data<3*rms]=0
            fluxList_3.append( np.sum(data)*0.2*omega  )


            data[data<4*rms]=0
            fluxList_4.append( np.sum(data)*0.2*omega  )


            data[data<5*rms]=0
            fluxList_5.append( np.sum(data)*0.2*omega  )

            data[data<6*rms]=0
            fluxList_6.append( np.sum(data)*0.2*omega  )

            data[data<7*rms]=0
            fluxList_7.append( np.sum(data)*0.2*omega  )






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

    def getRawCOFITS(self,calCode):
        """

        :param calCode:
        :return:
        """

        return self.dataPath+calCode+".fits"




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
            flux=np.sum(coValues)*0.2
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
            flux=np.sum(coValues)*0.2
            fluxList.append(flux)


        del cloudIndex
        del coValues


        return fluxList




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


    def getFillingFactorByCloudID(self, CODataRaw, labelSets,cleanDataList,calCode,ID, saveRow=None , drawFigure=False, useSigmaCut=False,printFluxList=False ):
        """
        :param ID:
        :return:
        """
        self.calCode=calCode
        #test an ID of 8570

        if useSigmaCut:
            fluxList = self.getFluxListByIDSigmaCut( labelSets,cleanDataList, ID  )

        else:
            fluxList = self.getFluxListByID(CODataRaw, labelSets,cleanDataList,calCode,ID)


        #to test memory
        #return

        if printFluxList:
            print "The flux list is "
            print fluxList

        wmsipFilling, cfaFilling,  fittingParaAndError=self.getFillingFactorAndDraw(self.smoothFactors*self.rawBeamSize,fluxList,calID=ID,drawFigure=drawFigure)
        #print wmsipFilling
        #save
        if saveRow!=None:
            #print "aaaaaaaaaaaaaaaaa"
            #tbIndexOfCloud= np.where(  saveTB["_idx"] ==ID   )

            para,paraError= fittingParaAndError
            saveRow[self.ffMWISPCol]=wmsipFilling
            saveRow[self.ffCfACol]=cfaFilling

            saveRow[self.aCol]=para[0]
            saveRow[self.bCol]=para[1]
            saveRow[self.cCol]=para[2]


            saveRow[self.aErrorCol]=paraError[0]
            saveRow[self.bErrorCol]=paraError[1]
            saveRow[self.cErrorCol]=paraError[2]

        return wmsipFilling, cfaFilling,  fittingParaAndError

    def getIndices(self,labelSets, choseID ):
        Z0, Y0, X0, values1D =labelSets
        cloudIndices = np.where(values1D == choseID)

        cX0 = X0[cloudIndices]
        cY0 = Y0[cloudIndices]
        cZ0 = Z0[cloudIndices]

        return tuple([cZ0, cY0, cX0])

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





    def calFFByID(self,calCode,ID,drawFigure=True, useSigmaCut=True  ):

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
        TBName = self.getRawBeamTBByCalcode(calCode)
        cleanTB = Table.read(TBName)

        ffTB=self.addFFColnames( cleanTB )


        if useSigmaCut:
            cleanDataList=self.getNoiseDataList()

        else:
            cleanDataList=self.getCleanDataList(calCode)
        rawCOFITS=self.getRawCOFITS(calCode)
        CODataRaw, COHeadRaw = doFITS.readFITS( rawCOFITS )

        cleanFITSRawBeam = self.getCleanFITSName(calCode, 1)
        cleanDataSM1,head=doFITS.readFITS(cleanFITSRawBeam)


        clusterIndex1D = np.where(cleanDataSM1 > 0)
        clusterValue1D = cleanDataSM1[clusterIndex1D]

        Z0, Y0, X0 = clusterIndex1D
        labelSets=[Z0, Y0, X0, clusterValue1D ]

        print self.getFillingFactorByCloudID(CODataRaw, labelSets, cleanDataList, calCode, ID, saveRow=None,  drawFigure=drawFigure, useSigmaCut=useSigmaCut, printFluxList=True)

    def getFluxColName(self,smFactor):

        return "fluxSM{:.1f}".format( smFactor )


    def getSmoothFluxCol(self,smFITS, TB,  labelSets, sigmaCut=2.6 ):
        """
        :return:
        """

        ####
        print "Extracting flux from ", smFITS
        dataSm,headSm= doFITS.readFITS(smFITS)

        smFactor = self.getSmoothFactor(smFITS)
        colName=self.getFluxColName(smFactor)

        TB[colName]=TB["peak"]*0


        widgets = ['Calculating filling factors: ', Percentage(), ' ', Bar(marker='0', left='[', right=']'), ' ', ETA(), ' ',
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
            fluxID=np.sum(coValues)*0.2

            eachRow[colName] = fluxID

            pbar.update(i)

        pbar.finish()


    def getFluxList(self,row):
        """

        :param row:
        :return:
        """

        fluxList=[]

        for eachSm in self.smoothFactors:
            colName=self.getFluxColName(eachSm)
            fluxList.append( row[colName])

        return np.asarray(fluxList)



    def calculateFillingFactor(self,TBFile,drawFigure=False):

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
            fluxList = self.getFluxList(eachRow)
            ID=eachRow[ self.idCol ]
            wmsipFilling, cfaFilling,  fittingParaAndError=self.getFillingFactorAndDraw(beamArray,fluxList,calID=ID,drawFigure=drawFigure)

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



    def getFFForEachCloud(self,calCode=None, drawFigure=False, useSigmaCut=True, calAllCloud=True ):

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
        else:
            calCode=self.calCode

        TBName = self.getSmoothAndNoiseFITSSingle(smFactor=1.0,noiseFactor=0.0,getCleanTBFile=True)


        cleanTB = Table.read(TBName)
        ffTB=self.addFFColnames( cleanTB )

        rawCOFITS=self.getRawCOFITS(calCode)
        CODataRaw, COHeadRaw = doFITS.readFITS( rawCOFITS )

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


        ffTB.write("fluxTB_"+os.path.basename(TBName) ,overwrite=True )

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
            cleanDataList=self.getCleanDataList(calCode)

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

            self.getFillingFactorByCloudID(CODataRaw, labelSets,cleanDataList,calCode, ID, saveRow=eachR, drawFigure=drawFigure, useSigmaCut=useSigmaCut)


        pbar.finish()
            #break
        if calAllCloud:
            ffTB.write( calCode+"FillingFactorTBAll.fit",overwrite=True )

        else:
            ffTB.write( calCode+"FillingFactorTB.fit",overwrite=True )






    def getFillingErrorAndError(self, beamSize , ffTB):

        """

        :param ffTB:
        :return:
        """
        x=  beamSize #arcmin
        a=ffTB[self.aCol]
        b=ffTB[self.bCol]
        c=ffTB[self.cCol]

        stdA = ffTB[self.aErrorCol ]
        stdB = ffTB[self.bErrorCol ]
        stdC = ffTB[self.cErrorCol ]

        fx=ffFunction(x,a,b,c)

        f0= ffFunction(0,a,b,c)
        varX = stdA**2* ( dya(x,a,b,c) )**2 + stdB**2* ( dyb(x,a,b,c) )**2 + stdC**2


        var0=  stdA**2* ( dya(0,a,b,c) )**2 + stdB**2* ( dyb(0,a,b,c) )**2 + stdC**2


        varFF=varX/f0**2 + (-fx/f0**2)**2*var0

        return  fx/f0  ,  np.sqrt( varFF )
        #beamSize


    def addMWISPFFerror(self, TB ):
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




    def drawFillingRelation(self,calCode,  fillingTB, drawCode="area"):
        """

        :param fillingTB:
        :return:
        """
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

            axFF.set_xlim([0, 20 ])
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


        axFF.set_ylim([-0.2, 1.2  ])


        axFF.set_ylabel("Filling factor")



        plt.savefig( "{}FFindividual_{}.png".format(calCode, drawCode ), bbox_inches='tight', dpi=600)


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
            drawX= TB["sum"]*0.2



        drawY=  TB[self.ffMWISPCol]

        yerr=  TB[self.ffMWISPErrorCol]

        if not showYError:
            yerr =  None
        import matplotlib.cm as cm
        import matplotlib as mpl

        Vs=TB["v_cen"]

        cmap = plt.cm.jet
        # norm = matplotlib.colors.BoundaryNorm(np.arange(0,30,0.1), cmap.N)
        normV = mpl.colors.Normalize(vmin=-6, vmax=30)
        m = plt.cm.ScalarMappable(norm=normV, cmap=cmap)
        v_color = np.array([(m.to_rgba(v)) for v in Vs])

        if color=="gray":
            #c =Vs, cmap=cmap,norm=normV,
            #ax.errorbar(drawX,   drawY , yerr= yerr ,    markersize=markerSize , linestyle='none', c= v_color ,     marker= markerType , capsize=0,  elinewidth=elinewidth , lw=1, label= label )
            #ax.errorbar(drawX,   drawY , yerr= yerr ,    markersize=markerSize , linestyle='none', c= Vs ,  cmap=cmap,norm=normV,   marker= markerType , capsize=0,  elinewidth=elinewidth , lw=1, label= label )
            sc = ax.scatter(drawX, drawY, c=Vs, cmap=cmap, norm=normV, s=3, facecolors='none', lw=0.5, marker="o")
            return sc



        else:

            ax.errorbar(drawX,   drawY , yerr= yerr ,    markersize=markerSize , linestyle='none', c= color ,    marker= markerType , capsize=0,  elinewidth=elinewidth , lw=1, label= label )


    def zzz(self):
        """

        :return:
        """

