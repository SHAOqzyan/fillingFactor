
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






class checkFillingFactor(object):

    rootPath="./"

    #saveFITSPath="/home/qzyan/WORK/diskMWISP/fillingFactorData/"
    saveFITSPath="./"

    dataPath= saveFITSPath+"data/"
    rawFITS=  dataPath + "G2650Local30.fits"

    #tmpPath= rootPath +"tmpFiles/"
    tmpPath= saveFITSPath +  "tmpFiles/"
    #########out
    figurePath=  rootPath +"figurePath/"



    codeOutCO12="OutCO12"
    outCO12FITS =dataPath+ "outCO12.fits"

    codeOutCO13="OutCO13"
    outCO13FITS =dataPath+ "outCO13.fits"

    codeOutCO18="OutCO18"
    outCO18FITS =dataPath+ "outCO18.fits"

    ############## Local
    codeLocalCO12="LocalCO12"
    localCO12FITS =dataPath+ "localCO12.fits"


    codeLocalCO13="LocalCO13"
    localCO13FITS =dataPath+ "localCO13.fits"


    codeLocalCO18="LocalCO18"
    localCO18FITS =dataPath+ "localCO18.fits"



    ########### Sgr

    codeSgrCO12= "sgrCO12"
    sgrCO12FITS =dataPath+ "sgrCO12.fits"

    codeSgrCO13 = "sgrCO13"
    sgrCO13FITS = dataPath+ "sgrCO13.fits"

    codeSgrCO18 = "sgrCO18"
    sgrCO18FITS = dataPath+ "sgrCO18.fits"

    ###### Scu

    codeScuCO12= "scuCO12"
    scuCO12FITS = dataPath+ "scuCO12.fits"

    codeScuCO13 = "scuCO13"
    scuCO13FITS = dataPath+ "scuCO13.fits"

    codeScuCO18 = "scuCO18"
    scuCO18FITS = dataPath+ "scuCO18.fits"


    outVrange=[-79,-6] #km/s
    localVrange=[-6,30] #km/s
    sgrVrange=[30,70] #km/s
    scuVrange=[70,120] #km/s



    #tmpPath=  "/home/qzyan/WORK/dataDisk/fillingFactorData/"

    #rawRMS = 0.5
    MWISPrmsCO12=0.49 #K

    MWISPrmsCO13=0.26  #K
    MWISPrmsCO18=0.28 #K


    #smoothFactors =  np.arange(1.0,15.5,0.5  ) #compared to beam size, this is also the distance factor, that push cloud away by this factor
    smoothFactors =  np.arange(1.0, 10.5, 0.5  ) #compared to beam size, this is also the distance factor, that push cloud away by this factor


    noiseFactors =     np.arange(1.0,10.5,0.5  )  # do not care about the noise Factores now


    #just for test
    #smoothFactors = np.arange(1.0, 2,   0.5)  # compared to beam size, this is also the distance factor, that push cloud away by this factor
    #noiseFactors = np.arange(1.0, 2, 0.5)  # compared to  raw data rms

    calCode=""
    noiseStr="noise"

    rawBeamSize =   52. / 60


    normNoise = mpl.colors.Normalize(vmin= min(noiseFactors) , vmax= max(noiseFactors) )
    noiseColor = plt.cm.ScalarMappable(norm=normNoise, cmap=plt.cm.jet)

    normSmooth = mpl.colors.Normalize(vmin= min(smoothFactors) * rawBeamSize , vmax= max(smoothFactors)*rawBeamSize )
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
    ffCfACol="fillingFactorCfa"

    aCol="para_a"
    bCol="para_b"
    cCol="para_c"




    aErrorCol="error_a"
    bErrorCol="error_b"
    cErrorCol="error_c"





    def __init__(self):
        pass



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


    def getTBFixNoiseFactor(self,noiseFactor=1.0, dbParameters="S2P4Con1" ):
        """

        :param noiseFactor:
        :return:
        """

        tbList=[]
        smList=[]
        for eachSF in self.smoothFactors:
            searchStr="{}*SmFactor*{}*noiseFactor*{}*{}*_Clean.fit".format(self.calCode,eachSF,noiseFactor,dbParameters)

            searchStr=os.path.join(self.tmpPath,searchStr)

            tbFile=glob.glob( searchStr )

            if len(tbFile) ==1:


                tbList.append(  Table.read( tbFile[0] ) )
                smList.append( eachSF )


            else:
                print "Multiple data found please check your data"
                print tbFile


        return np.asarray( tbList ) , np.asarray( smList )





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

        FITSrms = doFITS.getRMSFITS(FITSfile , ""  , returnRMSValue=True )

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

        rawBeamSize = self.rawBeamSize # 52. / 60,
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

    def getSmFITSFileList(self):

        fileList=[]


        for i in self.smoothFactors:

            searchStr =  "{}*SmFactor_{}.fits".format(self.calCode, i  )
            searchStr = os.path.join( self.tmpPath, searchStr)
            fitsName=glob.glob(searchStr)

            if len(fitsName )==1:
                fileList.append(fitsName[0] )
            else:
                print fitsName
                print "None of multiple fits found!"


        return  fileList




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



    def getAbsNoiseCleanTBLIst(self,absK=0.5, dbscanCode="dbscanS2P4Con1"):

        tbList= []
        tbFileList= []
        for eachSF in self.smoothFactors:
            fileName="{}*SmFactor_{}*{}absK{}_Clean.fit".format(self.calCode, eachSF ,absK, dbscanCode )

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
        wmsipFilling = ffFunction(self.rawBeamSize, params[0], params[1], params[2]) / ffFunction(0, params[0], params[1],   params[2])
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

        plt.savefig(self.figurePath+"{}.png".format( saveTag ), bbox_inches='tight', dpi=600)



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


    def getFluxListByID(self,CODataRaw, cleanDataSM1, cleanDataList,  calCode,ID  ):
        """

        :param ID:
        :return:
        """

        fluxList=[]

        #rawCOFITS=self.getRawCOFITS(calCode)

        #CODataRaw, COHeadRaw = doFITS.readFITS( rawCOFITS )

        #step 1, get the
        #data,head=doFITS.readFITS(cleanFITSRawBeam)
        cloudIndex= np.where( cleanDataSM1==ID )

        rawCOValues=  CODataRaw[cloudIndex]


        for  dataSm in cleanDataList :

            #print eachSm
            #cleanFITSSm =  self.getCleanFITSName(calCode, eachSm )

            #dataSm, headSm = doFITS.readFITS( cleanFITSSm )

            smLabels=  dataSm[cloudIndex]

            #indexSelect = np.logical_and( cloudMaskRaw , dataSm>0 )

            coValues=    rawCOValues[smLabels>0]
            flux=np.sum(coValues)*0.2

            fluxList.append(flux)

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


    def getFillingFactorByCloudID(self, CODataRaw, labelSets,cleanDataList,calCode,ID, saveTB=None , drawFigure=False, useSigmaCut=False ):
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



        wmsipFilling, cfaFilling,  fittingParaAndError=self.getFillingFactorAndDraw(self.smoothFactors*self.rawBeamSize,fluxList,calID=ID,drawFigure=drawFigure)

        #save
        if saveTB!=None:

            tbIndexOfCloud= np.where(  saveTB["_idx"] ==ID   )

            para,paraError= fittingParaAndError
            saveTB[tbIndexOfCloud][self.ffMWISPCol]=wmsipFilling
            saveTB[tbIndexOfCloud][self.ffCfACol]=cfaFilling

            saveTB[tbIndexOfCloud][self.aCol]=para[0]
            saveTB[tbIndexOfCloud][self.bCol]=para[1]
            saveTB[tbIndexOfCloud][self.cCol]=para[2]


            saveTB[tbIndexOfCloud][self.aErrorCol]=paraError[0]
            saveTB[tbIndexOfCloud][self.bErrorCol]=paraError[1]
            saveTB[tbIndexOfCloud][self.cErrorCol]=paraError[2]

    def getIndices(self,labelSets, choseID ):
        Z0, Y0, X0, values1D =labelSets
        cloudIndices = np.where(values1D == choseID)

        cX0 = X0[cloudIndices]
        cY0 = Y0[cloudIndices]
        cZ0 = Z0[cloudIndices]

        return tuple([cZ0, cY0, cX0])

    def getFFForEachCloud(self,calCode,drawFigure=False, useSigmaCut=False ):

        TBName = self.getRawBeamTBByCalcode(calCode)
        cleanTB = Table.read(TBName)

        ffTB=self.addFFColnames( cleanTB )
        self.calCode= calCode

        #CODataRaw, cleanDataSM1, cleanDataList,

        rawCOFITS=self.getRawCOFITS(calCode)
        CODataRaw, COHeadRaw = doFITS.readFITS( rawCOFITS )

        cleanFITSRawBeam = self.getCleanFITSName(calCode, 1)
        cleanDataSM1,head=doFITS.readFITS(cleanFITSRawBeam)

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
            gc.collect()

            ID=eachR["_idx"]
            pbar.update(i)

            self.getFillingFactorByCloudID(CODataRaw, labelSets,cleanDataList,calCode, ID,saveTB=ffTB,drawFigure=drawFigure, useSigmaCut=useSigmaCut)


        pbar.finish()
            #break
        ffTB.write( calCode+"FillingFactorTB.fit",overwrite=True )





    def zzz(self):
        """

        :return:
        """

