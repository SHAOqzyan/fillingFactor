#from fillingfactor import checkFillingFactor
from fillingfactor import *

import os
import sys
from astropy.table import Table,vstack
from myPYTHON import *
from matplotlib.colors import LogNorm

import fillingfactor as ff


doFITS=myFITS() #used to deal fits with myPYTHON

doFF=checkFillingFactor()

convertingRMS=0.49

def forward(x):
    return x*convertingRMS



def inverse(x):
    return  x/convertingRMS




def ffModelCutoff(x,a  ,b  ):


    processX=x-   b


    return   processX**2/( processX + a)**2  #processX /(processX**2 + a )  #1- np.exp(-b* processX )

    #return (x-2  )**2 / ((x-2 )**2 +a)
    #return  b*(np.arctan(a*(x  ) )  +1.5)

    #return  1- a*np.exp(-b*(x-2)**2)     #(  x- 2)**2/(x+a)**2
    #return   (x-2)**2/(x**2+a)    #(  x- 2)**2/(x+a)**2

class fillingMain(object):

    saveOverallCutoff= "/home/qzyan/WORK/myDownloads/fillingFactor/saveOverallCutoff/"

    ##############

    def __init__(self):
        pass


    def addNoiseToSmFiles(self,calCode,targetNoise):
        """

        :return:
        """

        doFF.calCode=calCode
        smFiles=doFF.getSmFITSFileList()

        for eachF in  smFiles:
            doFF.addNoiseSingle(eachF, absoluteNoise= targetNoise )


    def cleanByCode(self,calCode,targetNoise, removeFITS=True):
        """

        :param calCode:
        :param targetNoise:
        :return:
        """

        #get all noise fits
        doFF.calCode = calCode
        noiseFiles=doFF.getAbsNoiseFITSName( targetNoise )



        for eachF in noiseFiles:

            doFF.cleanFITSsigma2(eachF, removeFITS= removeFITS )


    def drawFillingFactor(self,calCode,drawFigure=False,inputID=None,testPoly=False,individuals=False):
        """
        :return:
        """

        doFF.calCode=calCode
        smTBFile = doFF.getSmoothAndNoiseFITSSingle(getCleanTBFile=True)



        fluxTBFile="fluxTB_"+os.path.basename(smTBFile)





        if not os.path.isfile( fluxTBFile):
            print "flux TB does not exist, please check your data or calcode"

            return

        print fluxTBFile

        #getFillingfactor
        ffTBName=doFF.calculateFillingFactor(fluxTBFile, drawFigure=drawFigure, inputID=inputID,testPoly=testPoly,individuals=individuals )




    def getFillingFactorAndEdge(self,calCode,drawFigure=False):

        #get fluxTB
        doFF.calCode=calCode
        smTBFile = doFF.getSmoothAndNoiseFITSSingle(getCleanTBFile=True)




        fluxTBFile="fluxTB_"+os.path.basename(smTBFile)

        if not os.path.isfile( fluxTBFile):
            print "flux TB does not exist, please check your data or calcode"

            return



        #getFillingfactor
        ffTBName=doFF.calculateFillingFactor(fluxTBFile, drawFigure=drawFigure)

        #get edgeInfo
        cleanFITSName = doFF.getSmoothAndNoiseFITSSingle(  getCleanFITS=True)

        doFF.getCloudCubes(cleanFITSName, ffTBName, writeFITS=False )


    def getFFTB(self,calCode):
        """
        return the table name of containning cloud Info
        :param calCode:
        :return:
        """
        doFF.calCode=calCode
        smTBFile = doFF.getSmoothAndNoiseFITSSingle(getCleanTBFile=True)

        return "pureEdgeInfo_fillingFactor_fluxTB_"+ os.path.basename( smTBFile )

    def getcutOFFFluxTB(self,calCode):

        doFF.calCode=calCode
        smTBFile = doFF.getSmoothAndNoiseFITSSingle(getCleanTBFile=True)

        return "cutoffFlux_pureEdgeInfo_fillingFactor_fluxTB_"+ os.path.basename( smTBFile )


    def cleanFITS(self,calCode,onlyFirst=False):

        doFF.calCode = calCode
        noiseFiles = doFF.getSmoothListFixNoise(noiseFactor=0.0)
        for eachF in noiseFiles:

            smFactor=doFF.getSmoothFactor(eachF)

            if onlyFirst:
                if smFactor>1.0:
                    continue

            doFF.cleanFITSsigma2(eachF)


    def removeUselessFITS(self,calCode):
        doFF.calCode = calCode
        for eachF in doFF.smoothFactors:
            if eachF==1.0:
                print "Kepp all files of smooth factor 1"
                continue
            print "Delete smooth factor {} for ".format(eachF),doFF.calCode
            rawDBSCANlabel=doFF.getSmoothAndNoiseFITSSingle(    smFactor= eachF  , getRawDBSCANFITS=True)
            cleanDBSCANlabel=doFF.getSmoothAndNoiseFITSSingle(   smFactor= eachF  , getCleanFITS=True)

            if rawDBSCANlabel is not None and os.path.isfile(rawDBSCANlabel):

                print "Removing ", rawDBSCANlabel
                os.remove(rawDBSCANlabel)

            if cleanDBSCANlabel is not None and os.path.isfile(cleanDBSCANlabel):
                print "Removing ", cleanDBSCANlabel
                os.remove(cleanDBSCANlabel)





    #def examineNoiseClouds(self,tbFile, ):

    def testBackChannelFFT(self,drawChannel=31):
        """

        :return:
        """

        testFITSFile = "/home/qzyan/WORK/diskMWISP/fillingFactorData/data/rawSgrCO12.fits"

        #drawChannel=31
        rawCOData,rawCOHead= myFITS.readFITS(testFITSFile)
        im=rawCOData[drawChannel]
        im = im.astype(float)
        im=np.nan_to_num(im)
        wcsCO = WCS(rawCOHead)

        #crop fits
        #im=im[:450,2200:]

        def plot_spectrum(im_fft):

            # A logarithmic colormap
            plt.imshow(np.abs(im_fft), norm=LogNorm(vmin=5))
            plt.colorbar()

        from scipy import fftpack
        im_fft = fftpack.fft2(im)

        # Show the results

        #def plot_spectrum(im_fft):
        from matplotlib.colors import LogNorm
        # A logarithmic colormap

            #plt.colorbar()
            #print "run"


        #print np.nanmax(im_fft)
        #print np.nanmin(im_fft)
        keep_fraction=0.2
        im_fft2 = im_fft.copy()

        # Set r and c to be the number of rows and columns of the array.
        r, c = im_fft2.shape

        # Set to zero all rows with indices between r*keep_fraction and
        # r*(1-keep_fraction):
        im_fft2[int(r * keep_fraction):int(r * (1 - keep_fraction))] = 0

        # Similarly with the columns:
        im_fft2[:, int(c * keep_fraction):int(c * (1 - keep_fraction))] = 0



        plt.figure()
        #plot_spectrum(im_fft2)
        im_new = fftpack.ifft2(im_fft2).real

        axRaw = pywcsgrid2.subplot(211, header=wcsCO)
        axProcess = pywcsgrid2.subplot(212, header=wcsCO)

        #im_new=im_new-im
        axRaw.imshow(im, plt.cm.gray,origin='lower',vmin=0.5,vmax=8)
        axProcess.imshow(im_new, plt.cm.gray,origin='lower',vmin=0.5,vmax=8)

        plt.title('Reconstructed Image')
        #plt.imshow(np.abs(im_fft), norm=LogNorm(vmin=10 ,vmax=6853))

        #plt.title('Fourier transform')
        plt.savefig(   "fftTest{}.png".format(drawChannel) , bbox_inches='tight', dpi=600)


    def getCombinCutoffTable(self):
        combineTB=None

        tableList=[ "rawLocalCO12BFFclean.fit" ,"rawLocalCO18BFFclean.fit","rawOutCO13BFFclean.fit","rawScuCO13BFFclean.fit",\
                    "rawSgrCO12BFFclean.fit","rawSgrCO18BFFclean.fit","rawLocalCO13BFFclean.fit","rawOutCO12BFFclean.fit",\
                    "rawScuCO12BFFclean.fit","rawScuCO18BFFclean.fit","rawSgrCO13BFFclean.fit"]


        for eachTB in tableList:

            testTB=Table.read(eachTB)

            if combineTB is None:
                    combineTB=testTB

            else:
                combineTB=vstack([combineTB,testTB])

        return combineTB


    def getCombinCutoffTableSurvey(self):
        combineTB=None

        tableList=[ "surveyCfACO12Q2BFFclean.fit" ,"surveyGRSCO13Q1BFFclean.fit","surveyOGSCO12Q2BFFclean.fit","surveyCOHRSCO32Q1BFFclean.fit" ]


        for eachTB in tableList:

            testTB=Table.read(eachTB)

            if combineTB is None:
                    combineTB=testTB

            else:
                combineTB=vstack([combineTB,testTB])

        return combineTB


    def getCombinFFTB(self,codeList):
        """

        :return:
        """
        combineTB=None

        for eachCode in  codeList:

        #for eachCode in doFF.allRawCodeList :

            doFF.calCode=  eachCode  #doFF.codeRawLocalCO18
            tb=doMain.getFFTB(doFF.calCode)

            testTB=Table.read(tb)
            print "Total number of clouds of ",eachCode, " is ",len(testTB)
            #add size in beam colName



            testTB[doFF.sizeInBeamCol] =doFF.getCloudSize(testTB)/doFF.getBeamSize()
            testTB[doFF.peakSigma] = doFF.getPeakSigma(testTB) #testTB["peak"] / doFF.getMeanRMS()

            if combineTB is None:
                    combineTB=testTB

            else:
                combineTB=vstack([combineTB,testTB])

        return combineTB



    def fittingFFWithAllTBs(self,codeList,saveTag="",drawCutOff=False):
        """
        fitting an overall filling factor, with all samples
        :return:
        """

        #according to the code list

        #for eachCoe in doFF.allRawCodeList:

        combineTB=self.getCombinFFTB(codeList)

        # use BeamFactorUnitFor combined clouds
        doFF.fittingAngularSizeFF( combineTB,showSizeRange=[-1,50],saveTag=saveTag,useBeamFactorUnit=True ,showColor="peak", drawCutOff = drawCutOff )



    def drawChiSquareTest(self,trainingRatio=0.8):
        """

        :return:
        """
        print "Testing chi-square of three functions..."

        chiSquareList1 = []
        chiSquareList2 = []
        chiSquareList3 = []



        errorCutList =  np.arange(0.1,2.1,0.1)  #[0.50,0.40,0.30,0.20,0.10,0.05]
        #errorCutList= [0.2, 0.5]

        doFF.calCode = doFF.codeRawLocalCO12

        if 1: #use single molecular cloud

            tb = self.getFFTB(doFF.calCode)
            ffTB = Table.read(tb)


            ffTB = doFF.addFFValues(ffTB)
            print len( ffTB)
            ###### #####test with all molecular clouds
        else:
            ffTB=self.getCombinFFTB(doFF.allRawCodeList)
            ffTB = doFF.addFFValues(ffTB)



        #used to draw a demonstrated figure with a maximum error of 50%
        trainHalfErrorCut = None
        testHalfErrorCut = None


        paraHalfErrorCutff1 = None
        paraErrorHalfErrorCutff1 = None

        paraHalfErrorCutff2 = None
        paraErrorHalfErrorCutff2 = None

        paraHalfErrorCutff3 = None
        paraErrorHalfErrorCutff3 = None



        testFunction1= ffModelConvolve
        testFunction2=   ffModelSquare
        testFunction3 =  testffAndSizeFunc2

        for eachEcut in  errorCutList:

            processTB=doFF.selectByErrorCut(ffTB,errorCut=eachEcut)

            trainingTB, testTB = doFF.splitTBIntoTrainingTest(processTB, trainingRatio=trainingRatio)



            if 1:

                cs1,para1,paraError1=doFF.testThreeFunctionsByTraining( testFunction1 ,   trainingTB,testTB, trainingRatio=trainingRatio, errorCut = eachEcut  )
                chiSquareList1.append(cs1)

                cs2,para2,paraError2=doFF.testThreeFunctionsByTraining( testFunction2 ,   trainingTB,testTB,  trainingRatio=trainingRatio,  errorCut = eachEcut  )
                chiSquareList2.append(cs2)


                cs3,para3,paraError3 =doFF.testThreeFunctionsByTraining( testFunction3 ,  trainingTB,testTB,  trainingRatio=trainingRatio,  errorCut = eachEcut  )
                chiSquareList3.append(cs3)




            if eachEcut==0.5:
                trainHalfErrorCut = trainingTB
                testHalfErrorCut = testTB
                ###############################
                paraHalfErrorCutff1 = para1
                paraErrorHalfErrorCutff1 = paraError1

                paraHalfErrorCutff2 = para2
                paraErrorHalfErrorCutff2 = paraError2

                paraHalfErrorCutff3 = para3
                paraErrorHalfErrorCutff3 = paraError3




        #draw three functions
        plt.clf()
        fig = plt.figure(figsize=(16, 7))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 17, 'serif': ['Helvetica']})
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]




        #draw a figure with
        axFF = plt.subplot(121 )

        print "Traing number of examples",len(trainHalfErrorCut)
        print "test number   of examples",len(testHalfErrorCut)


        drawXTrain = doFF.getCloudSize(trainHalfErrorCut)
        drawYTrain = trainHalfErrorCut[doFF.ffMWISPCol]
        yErrorTrain = trainHalfErrorCut[doFF.ffMWISPErrorCol]

        axFF.errorbar(drawXTrain,   drawYTrain ,  yerr= yErrorTrain ,    markersize=1 , fmt='o',  color='gray' , capsize=0.1,  elinewidth=0.5 , lw=1 ,alpha=0.8  ,label="Training data (50\% maximum relative BFF errors)",zorder=1 )


        drawXTest = doFF.getCloudSize(testHalfErrorCut)
        drawYTest = testHalfErrorCut[doFF.ffMWISPCol]
        yErrorTest = testHalfErrorCut[doFF.ffMWISPErrorCol]

        axFF.errorbar(drawXTest,   drawYTest ,  yerr= yErrorTest ,    markersize=1 , fmt='o',  color='black' , capsize=0.1,  elinewidth=0.5 , lw=1 ,alpha=0.8  ,label="Validation data (50\% maximum relative BFF errors)" ,zorder=2 )



        showSizeRange=[-1,100]
        axFF.set_ylim([-0.01, 1.1])
        axFF.set_xlim( showSizeRange )

        axFF.legend(loc= 4 ,fontsize=14)

        #draw four lines
        x = np.arange(-0.1, showSizeRange[1], 0.01)

        print paraHalfErrorCutff1,"blueline parameter"
        axFF.plot(x, testFunction1(x, *paraHalfErrorCutff1 ), "b-", lw=1 ,zorder=3)
        #axFF.plot(x, testFunction1(x,  40 ), "b-", lw=1 ,zorder=3)

        axFF.plot(x, testFunction2(x, *paraHalfErrorCutff2 ), "g--", lw=1 ,zorder=3)
        axFF.plot(x, testFunction3(x, *paraHalfErrorCutff3 ), "r-", lw=1 ,zorder=3,alpha=0.8)

        axFF.axvline(x=0, ls="--", color='black', lw=0.8)
        axFF.axhline(y=0,ls="--",color='black',lw=0.8)



        axTest = plt.subplot(122 )

        axTest.plot(  errorCutList, chiSquareList1,'.-',color='blue', markersize=7,lw=1,label=r"$\eta_{\mathrm{bf}} =\mathit l^2/\left(\mathit l^2+ a\right)$" )
        axTest.plot(  errorCutList, chiSquareList2,'.--',color='green', markersize=7,lw=1,label=r"$\eta_{\mathrm{bf}}= \mathit l^2/\left(\mathit l+ a\right)^2$" )

        #axTest.plot(  errorCutList, chiSquareList2,'.--',color='green', markersize=7,lw=1,label=r"$\mathit f=a\left(1-\exp\left(-b\mathit l\right)\right)$" )
        #axTest.plot(  errorCutList, chiSquareList3,'.-',color='red', markersize=7,lw=1,label=r"$\mathit f=  \frac{{2}}{{\pi}}\mathrm{{arctan}}\left(a \mathit l\right)$" ,alpha=0.8)
        axTest.plot(  errorCutList, chiSquareList3,'.-',color='red', markersize=7,lw=1,label=r"$\eta_{\mathrm{bf}} =a\left(1-\exp\left(-b\mathit l\right)\right)$" ,alpha=0.8)

        axTest.legend(loc=4)
        axTest.set_ylim(-0.02,0.12 )

        axTest.set_xlabel("Maximum relative error of filling factors")
        #axTest.set_ylabel("Chi-Square")
        axTest.set_ylabel("Weighted residual rms of test data")

        #print "{:.0f}\% test data".format(1-trainingRatio), "What the hell?????????"
        at = AnchoredText(r"{:.0f}\% test data".format((1-trainingRatio)*100), loc=2, frameon=False)
        axTest.add_artist(at)

        saveTag= "modelTestTestRatio{:.2f}".format(1-trainingRatio)

        axFF.set_xlabel("Angular size (arcmin)")
        #axTest.set_ylabel("Chi-Square")
        axFF.set_ylabel("The beam filling factor")


        plt.savefig(   saveTag+".png"  , bbox_inches='tight', dpi=600)
        plt.savefig(   saveTag+".pdf"  , bbox_inches='tight' )







    def getCutffFF(self,calCode, drawFigure= False ,cutByPeak=False ,dim=3 ):

        doFF.calCode=calCode

        rawCOFITS=doFF.getRawCOFITS()
        labelFITS= doFF.getSmoothAndNoiseFITSSingle( smFactor=1.0,  noiseFactor=0.0 , getCleanFITS=True)
        TBFile=self.getFFTB(calCode)
        rmsFITS=doFF.getRMSFITS()
        doFF.getSmoothFluxColCutoff(rawCOFITS,labelFITS,rmsFITS, TBFile ,drawFigure=drawFigure,cutByPeak= cutByPeak ,dim=dim   )







    def getFITSBySuffix(self,outPath,ID,suffix="PureInt.fits"):


        cloudName= "Cloud{}*{}".format(ID,suffix)
        searchStr= os.path.join( outPath,cloudName)

        return glob.glob(searchStr)[0]



    def getIntensityStd(self,TB):
        """

        :return:
        """
        intFITSPath="/home/qzyan/WORK/myDownloads/fillingFactor/intFITS/"

        stdList=[]
        for eachRow in TB:
            ID=eachRow["_idx"]

            intFITS=     self.getFITSBySuffix(intFITSPath, ID ,suffix="PureInt.fits")

            dataInt,headInt=doFITS.readFITS(intFITS)

            coValues = dataInt[dataInt>0]

            std=np.std(coValues,ddof=1)


            stdList.append( std )

        return  np.asarray(stdList)

    def getCloudRow(self,TB,id):
        """
        return the
        :param TB:
        :param id:
        :return:
        """

        idCol=TB["_idx"]

        return Table( TB[idCol==id] )[0]

    #def compareCloudSameSize(self,drawIDs=[380924,414408,81770, 258207 ]):
    def compareCloudSameSize(self,drawIDs=[43900,262017, 227912 ]):
        intFITSPath="/home/qzyan/WORK/myDownloads/fillingFactor/intFITS/"
        boxSize=0.8 #degree
        cloudTB= Table.read(  "testPeakFF.fit" )

        peakSigma = doFF.getPeakSigma(cloudTB)
        cloudTB[doFF.peakSigma] = peakSigma

        cloudSize=doFF.getCloudSize(cloudTB)

        cloudTB["size"] = cloudSize



        if 0:
            cloudTB.sort("peak")
            stdList=self.getIntensityStd(cloudTB)
            peakSigma=doFF.getPeakSigma(cloudTB)
            ffCol= cloudTB[doFF.ffMWISPCol]


            fig = plt.figure(1, figsize=(10 ,8))
            rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica'], "size": 15})
            # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

            rc('text', usetex=True)
            mpl.rcParams['text.latex.preamble'] = [
                r'\usepackage{tgheros}',  # helvetica font
                r'\usepackage{sansmath}',  # math-font matching  helvetica
                r'\sansmath'  # actually tell tex to use it!
                r'\usepackage{siunitx}',  # micro symbols
                r'\sisetup{detect-all}',  # force siunitx to use the fonts
            ]

            ax=plt.subplot(111)

            ax.scatter(  peakSigma , ffCol,s=10,color="blue"  )
            at = AnchoredText("$^{12}\mathrm{CO}$ clouds in the Local arm", loc=4, frameon=False, prop={"color": "black"},alpha=0.6)
            ax.add_artist(at)

            at = AnchoredText("$9\leq \mathit{l}\leq11$ arcmin", loc=1, frameon=False, prop={"color": "black"})
            ax.add_artist(at)



            ax.set_xlabel(r"Peak ($\sigma$)")
            ax.set_ylabel("The beam filling factor")
            plt.savefig(   "peakBFFelation.pdf", bbox_inches='tight')





        ####



        testID = drawIDs[2]
        #plt.clf()
        fig = plt.figure(1, figsize=(30 ,12))
        #fig = plt.figure( )

        rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica'], "size": 18})
        # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

        rc('text', usetex=True)
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        intFITS=self.getFITSBySuffix(intFITSPath,drawIDs[0],suffix="int.fits")


        dataCO,headCO=doFITS.readFITS(intFITS)


        WCScloud=WCS(headCO ,naxis=2)

        ##cloud 1
        testID = drawIDs[0]

        axCOCloud1 = pywcsgrid2.subplot(241, header=WCScloud)
        axBFFCloud1 = plt.subplot(245  )
        self.compareCloudSameSizeSingle( axCOCloud1 , axBFFCloud1, testID,cloudTB,boxSize,intFITSPath, fig)
        if 1:
            ##cloud 2
            testID = drawIDs[1]

            axCOCloud2 = pywcsgrid2.subplot(242, header=WCScloud)
            axBFFCloud2 = plt.subplot(246  )
            self.compareCloudSameSizeSingle( axCOCloud2 , axBFFCloud2, testID,cloudTB,boxSize,intFITSPath , fig)

            ##cloud 3
            testID = drawIDs[2]

            axCOCloud3 = pywcsgrid2.subplot(243, header=WCScloud)
            axBFFCloud3 = plt.subplot(247  )
            self.compareCloudSameSizeSingle( axCOCloud3 , axBFFCloud3, testID,cloudTB,boxSize,intFITSPath, fig)


            ##cloud 4
            #testID = drawIDs[3]

            #axCOCloud3 = pywcsgrid2.subplot(244, header=WCScloud)
            #axBFFCloud3 = plt.subplot(248  )
            #self.compareCloudSameSizeSingle( axCOCloud3 , axBFFCloud3, testID,cloudTB,boxSize,intFITSPath, fig)




        saveTag="cloudSameSize"
        plt.savefig(   saveTag+".png"  , bbox_inches='tight', dpi=600)
        plt.savefig(   saveTag+".pdf"  , bbox_inches='tight' )


    def compareCloudSameSizeSingle(self, axCOCloud1 , axBFFCloud1, testID,cloudTB,boxSize,intFITSPath, fig):

        #intFITS=self.getFITSBySuffix(intFITSPath, testID ,suffix="PureInt.fits")



        intFITS=self.getFITSBySuffix(intFITSPath, testID ,suffix="int.fits")
        maskFITS=self.getFITSBySuffix(intFITSPath, testID ,suffix="mask.fits")

        dataCO,headCO=doFITS.readFITS(intFITS)

        dataMask,headMask =doFITS.readFITS(maskFITS)

        WCScloud=WCS(headCO ,naxis=2)


        cloudRow=self.getCloudRow(cloudTB,testID)
        centerL=cloudRow["x_cen"]
        centerV= cloudRow["v_cen"]
        centerB=cloudRow["y_cen"]

        lRange=[ centerL-boxSize/2. , centerL+boxSize/2. ]
        bRange=[ centerB-boxSize/2. , centerB+boxSize/2. ]
        #img=axCOCloud1.imshow(dataCO, plt.cm.bone,origin='lower',vmin=0,vmax=10)

        cmapCO = plt.cm.bone
        cmapCO.set_bad('black')
        img=axCOCloud1.imshow(dataCO, cmap=cmapCO,origin='lower',norm = LogNorm(vmin=1, vmax=20) )








        axCOCloud1[ WCScloud ].contour(dataMask, [1], colors='red',  linewidths=0.8)

        axCOCloud1.set_ticklabel_type("absdeg", "absdeg")
        axCOCloud1.axis[:].major_ticks.set_color("w")
        axCOCloud1.set_xlabel(r"Galactic Longitude ($^{\circ}$)")
        axCOCloud1.set_ylabel(r"Galactic Latitude ($^{\circ}$)")


        cloudName= doFITS.getCloudNameByLB(centerL,centerB)

        at = AnchoredText(cloudName  , loc=1, frameon=True, prop={"color":"black" } )
        axCOCloud1.add_artist(at )
        at = AnchoredText( "{:.1f} $\mathrm{{km}}\ \mathrm{{s}}^{{-1}}$".format(centerV), loc=3 , frameon=False, prop={"color":"white"} )
        axCOCloud1.add_artist(at)

        leftBottomCorner = WCScloud.wcs_world2pix( max(lRange),min(bRange),0 )
        rightTopCorner = WCScloud.wcs_world2pix( min(lRange),max(bRange),0 )

        if 0: #colorbar
            from matplotlib.axes import Axes
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(axCOCloud1)
            cax1 = divider.new_horizontal("5%", pad=0.05, axes_class=Axes)
            fig.add_axes(cax1)
            #cax1 = divider.append_axes("right", size="3%", pad=0.05)
            cb = fig.colorbar(img, cax=cax1)
            cax1.set_ylabel(r"Intensity ($\rm K\  km\ s^{-1}$)")

        if 1: #colorbar
            from mpl_toolkits.axes_grid1.inset_locator import inset_axes

            from matplotlib.axes import Axes
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            #divider = make_axes_locatable(axCOCloud1)
            #cax1 = divider.new_horizontal("5%", pad=0.05, axes_class=Axes)
            cax1 = inset_axes(axCOCloud1,
                             width="3%",  # width = 10% of parent_bbox width
                             height="100%",  # height : 50%
                             loc=3,
                             bbox_to_anchor=(1.01, 0, 1, 1),
                             bbox_transform=axCOCloud1.transAxes,
                             borderpad=0.
                             )
            fig.add_axes(cax1)

            cb = fig.colorbar(img, cax=cax1)
            cax1.set_ylabel(r"Intensity ($\rm K\  km\ s^{-1}$)")

            tickesArray = np.asarray([1, 2,4,8, 16 ])
            # tickesArray=tickesArray**2
            cb.set_ticks([ 1, 2,4,8, 16 ])
            #cb.set_ticklabels(["A", "B", "C", "D" , "D"])
            #cb.ax.set_yticks( tickesArray )
            #cb.ax.set_yticklabels(map(str, tickesArray))
            cb.set_ticklabels( map(str, tickesArray)  )

        ############### beam filling factor
        fluxList,fluxError=doFF.getFluxList(cloudRow)
        processBeam = doFF.getBeamSize() #in arcmin

        beamArray= doFF.smoothFactors* processBeam

        axBFFCloud1.errorbar(beamArray,   fluxList , yerr= fluxError ,    markersize=2 , fmt='o', c= 'red' ,     capsize=0.5,  elinewidth=1 , lw=1 )

        a,b,c=cloudRow[doFF.aCol], cloudRow[doFF.bCol], cloudRow[doFF.cCol]
        fittingX = np.arange(0, np.max(beamArray), 0.01)

        bff=cloudRow[ doFF.ffMWISPCol ]
        pcov=doFF.getPcov(cloudRow)

        bff,bffError=doFF.getFFandError( (a,b,c),pcov ,doFF.getBeamSize() )


        axBFFCloud1.plot(fittingX, ffFunction(fittingX, a,b,c), color='blue', lw=1.0 ,label=  " $\eta_{{\mathrm{{bf}}}} = {:.3f}\pm{:.3f}$".format( bff, bffError ) )
        at = AnchoredText( cloudName  , loc=3, frameon=False, prop={"color":"black"} )
        axBFFCloud1.add_artist(at)


        sizeInt=int(cloudRow["size"]  )

        sizeF1=  int(  (cloudRow["size"] -sizeInt )*10   )

        at = AnchoredText( r"{:.1f}$\sigma$ mean SNR ($\mathit l$ =  {}$.\mkern-4mu^\prime${})".format(cloudRow["meanSNR"],sizeInt, sizeF1)  , loc=5, frameon=False, prop={"color":"black"} )
        axBFFCloud1.add_artist(at)



        axBFFCloud1.set_xlabel("Beam size (arcmin)")
        #axFitting.set_ylabel("Flux (arcmin)")
        axBFFCloud1.set_ylabel(r"Flux ($\rm K\ km\ s$$^{-1}$ $\mathrm{\Omega}_\mathrm{A}$)")
        axBFFCloud1.legend(loc=1, handlelength = 0.5 )



        axBFFCloud1.axvline(x=0, ls="--", color='black',lw=1)

        axBFFCloud1.set_xlim(-0.2, 9  )

        #plt.subplots_adjust(wspace=0.1)
        plt.subplots_adjust(wspace=0.3)
        axCOCloud1.set_xlim( leftBottomCorner[0], rightTopCorner[0] )
        axCOCloud1.set_ylim( leftBottomCorner[1], rightTopCorner[1] )


    def testPeakBFF(self):
        """
        examine the beam filling factor calculated with peak values
        :return:
        """
        doFF.calCode =  doFF.codeRawLocalCO13
        tbFile = "fillingFactor_fluxTB_rawLocalCO13rawLocalCO13_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"  #doMain.getFFTB(doFF.codeRawLocalCO13)
        FFTB = Table.read(tbFile)

        bffPeak,bffPeakError=doFF.calculateFillingFactorPeak(FFTB,drawFigure=False,inputID=None)

        FFTB[doFF.ffPeakMWISPCol] = bffPeak
        FFTB[doFF.ffPeakMWISPErrorCol] = bffPeakError

        ###draw
        beamArray=doFF.smoothFactors*doFF.getBeamSize()
        fig = plt.figure(1, figsize=(10, 8))
        rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica'], "size": 15})
        # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

        rc('text', usetex=True)
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        ax = plt.subplot(111)

        cloudSize=doFF.getCloudSize(FFTB )
        ax.scatter(  FFTB[doFF.ffMWISPCol] ,FFTB[doFF.ffPeakMWISPCol] ,color='blue' ,s=2  )

        ax.plot( [0,1],[0,1], '--',color= "red" ,lw=1 )


        ax.set_xlabel(r"Beam filling factor with total flux")
        ax.set_ylabel("Beam filling factor with Peak")
        plt.savefig("testPeakBFF.pdf", bbox_inches='tight')
        plt.savefig("testPeakBFF.png", bbox_inches='tight',dpi=600)


    def produceSameSizeClouds(self):
        """

        :return:
        """


        outPath = "/home/qzyan/WORK/myDownloads/fillingFactor/intFITS"

        # from mwispDBSCAN import MWISPDBSCAN
        doDBSCAN = MWISPDBSCAN()

        doFF.calCode = doFF.codeRawLocalCO12

        tbFile = doMain.getFFTB(doFF.codeRawLocalCO12)
        tb = Table.read(tbFile)

        size = doFF.getCloudSize(tb)

        selectionCriteria = np.logical_and(size >= 9, size <= 11)  # arcmin

        newTB = tb[selectionCriteria]

        newTB = newTB[newTB[doFF.ffMWISPCol] > 0]
        newTB = newTB[newTB[doFF.ffMWISPCol] < 1]

        print len(newTB), "total number of test clouds"

        newTB["meanSNR"] = newTB['sum']/newTB['pixN']/doFF.getMeanRMS()



        newTB.write("testPeakFF.fit", overwrite=True)

        return

        doDBSCAN.produceCloudIntFITS(doFF.getRawCOFITS(),
                                     "/home/qzyan/WORK/diskMWISP/fillingFactorData/tmpFiles/rawLocalCO12rawLocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fits",
                                     newTB, \
                                     outputPath=outPath,   useTB=True)



    def drawCutOFFBFFSurvey(self,calCode ="survey",saveTag="cutoffBFFSurvey",showSizeRange=[2,10] ,inputTB=None,errorThresh=0.2, showLegend=False, cloudSource="" ):
        """
        draw survey code
        :return:
        """



        drawTB=self.getCombinCutoffTableSurvey()
        inputTB=drawTB
        # BFFF,BffERROR= doFF.getFillingAndError(  doFF.getBeamSize() , drawTB, useCutoff=False  )

        # remove bad data

        print len(drawTB)
        drawTB = drawTB[drawTB["cutFFfillingFactorMWISP"] > 0]
        drawTB = drawTB[drawTB["cutFFfillingFactorMWISP"] <= 1]

        drawTB = drawTB[drawTB["cutFFfillingFactorMWISPError"] > 0]
        drawTB = drawTB[drawTB["cutFFpara_a"] > 0]
        drawTB = drawTB[drawTB["cutFFpara_sigma"] > 0]
        drawTB = drawTB[drawTB["cutFFpara_mu"] > 0]
        drawTB = drawTB[drawTB["cutFFpara_b"] > 0]
        drawTB = drawTB[drawTB["cutFFpara_c"] > 0]

        print len(drawTB), "after"

        if errorThresh is not None:
            drawTBGood = drawTB[
                drawTB["cutFFfillingFactorMWISPError"] / drawTB["cutFFfillingFactorMWISP"] < errorThresh]
            print len(drawTBGood), "after error threshold"

        fig = plt.figure(figsize=(10, 8))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 22, 'serif': ['Helvetica']})

        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        ffTB = doFF.cleanCutoffTB(drawTB)

        axFFvel = fig.add_subplot(1, 1, 1)

        if 1:
            drawX = ffTB[doFF.meanSNRcol]  # ffTB["sum"]/ffTB["pixN"]   # ffTB["peak"]  # self.getCloudSize(useTB)
            # drawX =   doFF.getCloudSize( ffTB  ) #ffTB["peak"]  # self.getCloudSize(useTB)
            # drawX = ffTB["pixN"]

            # drawX =  ffTB["sum"] # ffTB["sum"]/ffTB["pixN"]   # ffTB["peak"]  # self.getCloudSize(useTB)

            drawY = ffTB[doFF.cutFFffMWISPCol]
            yError = ffTB[doFF.cutFFffMWISPErrorCol]

            para, paraError = doFF.fittingFFAndSize(ffModelCutoff, drawX, drawY, yError, sizeError=None, useODR=False)
            # a,b=para
            print para, "(a,b) parameters and errors below", calCode
            print paraError
        else:

            drawX = doFF.getCloudSize(ffTB)  # ffTB["peak"]  # self.getCloudSize(useTB)
            doFF.addMWISPFFerror(ffTB)

            drawY = ffTB[doFF.ffMWISPCol] / ffTB[doFF.cutFFffMWISPCol]
            yError = ffTB[doFF.ffMWISPErrorCol]

        # only draw good clouds

        drawX = drawTBGood[doFF.meanSNRcol]  # ffTB["sum"]/ffTB["pixN"]   # ffTB["peak"]  # self.getCloudSize(useTB)
        drawY = drawTBGood[doFF.cutFFffMWISPCol]
        yError = drawTBGood[doFF.cutFFffMWISPErrorCol]

        ffTB = drawTBGood

        # curve fitting
        if inputTB is None:
            ffTB.write(calCode + "BFFclean.fit")
        ####get color bar
        Vs = ffTB["v_cen"]
        normV = mpl.colors.Normalize(vmin=np.min(Vs), vmax=np.max(Vs))
        cmap = plt.cm.jet
        m = plt.cm.ScalarMappable(norm=normV, cmap=cmap)
        #axFFvel.scatter(drawX, drawY, edgecolors="blue", s=2, facecolors='none', lw=0.5, marker="o", label="", zorder=0,  alpha=0.8)
        # aa,bb,cc=axFFvel.errorbar(drawX, drawY, yerr=yError, markersize=0.0, fmt='o',  color="blue" ,   capsize=0.0,   elinewidth=0.7, lw=0.5, alpha=0.8, label="",  zorder=1)
        # cc[0].set_color(  m.to_rgba(Vs) )



        armStr = doFF.getArmStr()
        lineStr = doFF.getLineStr()
        fillingWMISP = r"{} molecular clouds in the {} arm".format(lineStr, armStr)


        #at = AnchoredText("", loc=4, frameon=False)
        #axFFvel.add_artist(at)
        #tableList=[ "surveyCfACO12Q2BFFclean.fit" ,"surveyGRSCO13Q1BFFclean.fit","surveyOGSCO12Q2BFFclean.fit","surveyCOHRSCO32Q1BFFclean.fit" ]


        ################# draw  COHRS
        if 1:

            cfaFFTBFile="surveyCOHRSCO32Q1BFFclean.fit" ###
            cfaFFTB= Table.read(cfaFFTBFile) ####
            cfaFFTB = cfaFFTB[  cfaFFTB["cutFFfillingFactorMWISPError"] / cfaFFTB["cutFFfillingFactorMWISP"] < errorThresh]

            drawX = cfaFFTB[doFF.meanSNRcol]  # ffTB["sum"]/ffTB["pixN"]   # ffTB["peak"]  # self.getCloudSize(useTB)
            drawY = cfaFFTB[doFF.cutFFffMWISPCol]

            axFFvel.scatter(drawX, drawY, edgecolors="blue", s=2.2, facecolors='none', lw=0.6, marker="o",     alpha=0.8, label=r"COHRS $^{12}\mathrm{CO}\left(J=3\rightarrow2\right)$",zorder=0 )



        ################# draw  OGSCO12
        if 1:

            cfaFFTBFile="surveyOGSCO12Q2BFFclean.fit" ###
            cfaFFTB= Table.read(cfaFFTBFile) ####
            cfaFFTB = cfaFFTB[  cfaFFTB["cutFFfillingFactorMWISPError"] / cfaFFTB["cutFFfillingFactorMWISP"] < errorThresh]

            drawX = cfaFFTB[doFF.meanSNRcol]  # ffTB["sum"]/ffTB["pixN"]   # ffTB["peak"]  # self.getCloudSize(useTB)
            drawY = cfaFFTB[doFF.cutFFffMWISPCol]

            axFFvel.scatter(drawX, drawY, edgecolors="red", s=2.2, facecolors='none', lw=0.6, marker="o",     alpha=0.8,label=r"OGS $^{12}\mathrm{CO}\left(J=1\rightarrow0\right)$"  ,zorder=1)



        ################# draw  GRS
        if 1:

            cfaFFTBFile="surveyGRSCO13Q1BFFclean.fit" ###
            cfaFFTB= Table.read(cfaFFTBFile) ####
            cfaFFTB = cfaFFTB[  cfaFFTB["cutFFfillingFactorMWISPError"] / cfaFFTB["cutFFfillingFactorMWISP"] < errorThresh]

            drawX = cfaFFTB[doFF.meanSNRcol]  # ffTB["sum"]/ffTB["pixN"]   # ffTB["peak"]  # self.getCloudSize(useTB)
            drawY = cfaFFTB[doFF.cutFFffMWISPCol]

            axFFvel.scatter(drawX, drawY, edgecolors="purple", s=2.2, facecolors='none', lw=0.6, marker="o",     alpha=0.8,label=r"GRS $^{13}\mathrm{CO}\left(J=1\rightarrow0\right)$" ,zorder=2)



        ################# draw CfA
        if 1:

            cfaFFTBFile="surveyCfACO12Q2BFFclean.fit" ###
            cfaFFTB= Table.read(cfaFFTBFile) ####
            cfaFFTB = cfaFFTB[  cfaFFTB["cutFFfillingFactorMWISPError"] / cfaFFTB["cutFFfillingFactorMWISP"] < errorThresh]

            drawX = cfaFFTB[doFF.meanSNRcol]  # ffTB["sum"]/ffTB["pixN"]   # ffTB["peak"]  # self.getCloudSize(useTB)
            drawY = cfaFFTB[doFF.cutFFffMWISPCol]

            axFFvel.scatter(drawX, drawY, edgecolors="green", s=2.2, facecolors='none', lw=0.6 , marker="o",     alpha=0.8,label=r"CfA 1.2-m  $^{12}\mathrm{CO}\left(J=1\rightarrow0\right)$" ,zorder=3)







        axFFvel.set_xlim(showSizeRange)
        axFFvel.set_ylim(0, 1)

        axFFvel.set_ylabel(r"$\mathit f$")  # ("The beam filling factor")
        # axFF.set_xlabel("Angular size (arcmin)")
        # axFFvel.set_xlabel("Angular size (arcmin)")
        #axFFvel.set_xlabel(r"Mean SNR ($\sigma$)")
        #axFFvel.set_xlabel(r"Mean SNR ($\sigma$)")
        axFFvel.set_xlabel(r"Mean SNR ($\sigma$)")

        

        ########
        para= [0.423, 2.329 ]
        # x = np.arange( 2.1875, showSizeRange[1], 0.01)
        if len(para) == 2:
            x = np.linspace(para[-1], showSizeRange[1], 1000)

        else:
            x = np.linspace(2.1875, showSizeRange[1], 1000)

        a, b = para
        formulaTex = r"$\mathit f = \frac{{ \left( \mathit x - {:.3f} \right) ^2}}{{\left( \left( \mathit x  - {:.3f} \right)+ {:.3f}  \right)^2 }}$".format(
            b, b, a)

        axFFvel.plot(x, ffModelCutoff(x, *para), "-", color='black', lw=1, zorder=2, label=formulaTex)
        ##########
        #if showLegend:
        axFFvel.legend(loc=5 , handlelength= 0.6,fontsize=18)

        plt.savefig(doFF.paperFigurePath + "{}FittingSNRAndBFF{}.png".format(calCode, saveTag), bbox_inches='tight',
                    dpi=300)
        plt.savefig(doFF.paperFigurePath + "{}FittingSNRAndBFF{}.pdf".format(calCode, saveTag), bbox_inches='tight')

    ##############

    def drawCutoffBFF(self,calCode ,saveTag="",showSizeRange=[2,10] ,inputTB=None,errorThresh=None, showLegend=True, cloudSource=None ,dim=5 ,showX2=False):
        """
        draw beam filling factor  caused by cutoff
        :return:
        """

        ffTB=self.getFFTB(calCode)

        doFF.calCode = calCode
        if inputTB is None:
            drawTB=Table.read("cutoffFF_{}".format(dim)+  ffTB ) #pureEdgeInfo_fillingFactor_fluxTB_rawLocalCO12rawLocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit")
        else:
            drawTB =  inputTB #pureEdgeInfo_fillingFactor_fluxTB_rawLocalCO12rawLocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit")


        #BFFF,BffERROR= doFF.getFillingAndError(  doFF.getBeamSize() , drawTB, useCutoff=False  )

        #remove bad data


        print len(drawTB),"before?"
        drawTB= drawTB [ drawTB["cutFFfillingFactorMWISP"] >0 ]

        drawTB= drawTB [ drawTB["cutFFfillingFactorMWISP"] <=1 ]

        drawTB= drawTB [ drawTB["cutFFfillingFactorMWISPError"] >0 ]

        #if dim==3:

        drawTB= drawTB [ drawTB["cutFFpara_a"] >0 ]

        if dim==5:

            drawTB= drawTB [ drawTB["cutFFpara_sigma"] >0 ]
            #drawTB= drawTB [ drawTB["cutFFpara_mu"] >0 ]
            drawTB= drawTB [ drawTB["cutFFpara_b"] >0 ]
            #drawTB= drawTB [ drawTB["cutFFpara_c"] >0 ]

        #what about relative error lesst than 20%
 

        



        print len(drawTB),"after"


        if errorThresh is not None:

            drawTBGood = drawTB[drawTB["cutFFfillingFactorMWISPError"] /drawTB["cutFFfillingFactorMWISP"] <errorThresh ]
            drawTB=drawTBGood
            print len(drawTBGood),"after error threshold"

        fig = plt.figure(figsize=(10, 8)  )
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size':  22 , 'serif': ['Helvetica']})

        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        ffTB = doFF.cleanCutoffTB( drawTB )


        axFFvel =  fig.add_subplot(1, 1, 1 )
        ax2 = axFFvel.twiny()
        if 1:
            drawX =  ffTB[doFF.meanSNRcol] # ffTB["sum"]/ffTB["pixN"]   # ffTB["peak"]  # self.getCloudSize(useTB)
            #drawX =   doFF.getCloudSize( ffTB  ) #ffTB["peak"]  # self.getCloudSize(useTB)
            #drawX = ffTB["pixN"]

            #drawX =  ffTB["sum"] # ffTB["sum"]/ffTB["pixN"]   # ffTB["peak"]  # self.getCloudSize(useTB)

            drawY = ffTB[doFF.cutFFffMWISPCol]
            yError = ffTB[doFF.cutFFffMWISPErrorCol]




            para,paraError=doFF.fittingFFAndSize( ffModelCutoff , drawX,drawY, yError,sizeError=None,  useODR=False )
            #a,b=para
            print para ,"(a,b) parameters and errors  of ", calCode
            print paraError
        else:

            drawX =   doFF.getCloudSize( ffTB  ) #ffTB["peak"]  # self.getCloudSize(useTB)
            doFF.addMWISPFFerror(ffTB )

            drawY = ffTB[doFF.ffMWISPCol]/  ffTB[doFF.cutFFffMWISPCol]
            yError = ffTB[doFF.ffMWISPErrorCol]


        #only draw good clouds
        if errorThresh is not None:
            drawX = drawTBGood[doFF.meanSNRcol]  # ffTB["sum"]/ffTB["pixN"]   # ffTB["peak"]  # self.getCloudSize(useTB)
            drawY = drawTBGood[doFF.cutFFffMWISPCol]
            yError = drawTBGood[doFF.cutFFffMWISPErrorCol]

            ffTB=drawTBGood

        #curve fitting
        if inputTB is None:
            ffTB.write(calCode+"BFFclean.fit",overwrite=True)
        ####get color bar
        Vs=ffTB["v_cen"]
        normV = mpl.colors.Normalize(vmin=np.min(Vs), vmax=np.max(Vs) )
        cmap = plt.cm.jet
        m = plt.cm.ScalarMappable(norm=normV, cmap=cmap)
        axFFvel.scatter(drawX, drawY, edgecolors =  "blue",  s= 2,    facecolors='none',  lw=0.5, marker="o", label="",zorder=0,alpha=0.8)
        #aa,bb,cc=axFFvel.errorbar(drawX, drawY, yerr=yError, markersize=0.0, fmt='o',  color="blue" ,   capsize=0.0,   elinewidth=0.7, lw=0.5, alpha=0.8, label="",  zorder=1)
        #cc[0].set_color(  m.to_rgba(Vs) )





        axFFvel.set_xlim(  showSizeRange  )
        axFFvel.set_ylim( 0,1 )

        axFFvel.set_ylabel(r"$\eta_{{\rm sns}}$") #("The beam filling factor")
        #axFF.set_xlabel("Angular size (arcmin)")
        #axFFvel.set_xlabel("Angular size (arcmin)")
        axFFvel.set_xlabel(r"Mean SNR ($\sigma$)")

        armStr=doFF.getArmStr()
        lineStr=doFF.getLineStr()
        fillingWMISP = r"{} molecular clouds in the {} arm".format(lineStr, armStr )


        if cloudSource is None:

            at = AnchoredText(fillingWMISP, loc=4, frameon=False)
            axFFvel.add_artist(at)
            
        else:
            at = AnchoredText(cloudSource , loc=4, frameon=False)
            axFFvel.add_artist(at)


        ########

        #x = np.arange( 2.1875, showSizeRange[1], 0.01)
        if len(para)==2:
            x = np.linspace(para[-1], showSizeRange[1],1000)
            a,b=para
        else:
            x = np.linspace( 2.1875, showSizeRange[1],1000)
            a =para[0]

            b=2.1875
        formulaTex = r"$\eta_{{\rm sns}} = \frac{{ \left( \mathit x - {:.3f} \right) ^2}}{{\left( \left( \mathit x  - {:.3f} \right)+ {:.3f}  \right)^2 }}$".format(b,b,a   )

        axFFvel.plot(x,  ffModelCutoff(x, *para), "-", color='black', lw=1,   zorder=2,label=formulaTex)
        ##########
        if showLegend:
            axFFvel.legend(loc=5 , handlelength=1.0)


        if showX2:
            ax2.set_xlim( showSizeRange )
            rangeofTemperature = np.asarray(showSizeRange)*doFF.getMeanRMS()

            minTem = int( min(rangeofTemperature) )   -1
            maxTem = int( max(rangeofTemperature) )+2
            dT = 0.2

            if "12" in calCode:
                dT=0.5

            tmpList= np.arange( minTem, maxTem+ dT, dT  )
            showSNR= tmpList /doFF.getMeanRMS()
            new_tick_locations= showSNR

            #get ticklabes


            print new_tick_locations,"?????????????????"

            ax2.set_xticks(new_tick_locations)
            ax2.set_xticklabels( self.getStr(tmpList)  )
            ax2.set_xlabel(r"Mean brightness temperature (K)")
            #secax = axFFvel.secondary_xaxis('top', functions=(forward, inverse ) )
            #secax.set_xlabel('Kelvin (K)')

            ax2.set_xlim( showSizeRange )
        plt.savefig(doFF.paperFigurePath+"{}FittingSNRAndBFF{}_dim{}.png".format(calCode,saveTag,dim), bbox_inches='tight', dpi=300)
        plt.savefig(doFF.paperFigurePath+"{}FittingSNRAndBFF{}_dim{}.pdf".format(calCode,saveTag,dim), bbox_inches='tight'  )



        #print   drawTB[doFF.cov00Col]

    def getStr(self,array):
        returnList=[]

        for eachDigiti in array:

            returnList.append( "{:.1f}".format( eachDigiti ) )

        return returnList

    def drawBFFrelatin(self,drawTB):
        """

        :param drawTB:
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

        #ffTB = doFF.cleanCutoffTB( drawTB )

        doFF.addMWISPFFerror(drawTB)
        drawTB=drawTB[ drawTB["fillingFactorMWISP"]>0  ]
        drawTB=drawTB[ drawTB["fillingFactorMWISP"]<1   ]
        drawTB=drawTB[ drawTB["fillingFactorMWISPError"] / drawTB["fillingFactorMWISP"]  <0.2   ]
        drawTB=drawTB[ drawTB["cutFFfillingFactorMWISPError"] / drawTB["cutFFfillingFactorMWISP"]  <0.2   ]


        #print drawTB.colnames
        axFFvel =  fig.add_subplot(1, 1, 1 )

        eta=drawTB["fillingFactorMWISP"]
        etaError=drawTB["fillingFactorMWISPError"]



        f=drawTB["cutFFfillingFactorMWISP"]
        fError=drawTB["cutFFfillingFactorMWISPError"]

        #aa,bb,cc=axFFvel.errorbar(eta, f, yerr=fError,xerr=etaError, markersize=0.5 , fmt='o',  color="blue" ,   capsize=0.0,   elinewidth=0.7, lw=0.5, alpha=0.8, label="",  zorder=1)
        aa,bb,cc=axFFvel.errorbar(  f,eta, yerr=etaError ,xerr=fError, markersize=0.5 , fmt='o',  color="blue" ,   capsize=0.0,   elinewidth=0.7, lw=0.5, alpha=0.8, label="",  zorder=1)
        #axFFvel.scatter(f, eta, edgecolors =  "blue",  s= 2,    facecolors='none',  lw=0.5, marker="o", label="",zorder=0,alpha=0.8)

        axFFvel.set_xlabel(r"$\mathit f$")
        axFFvel.set_ylabel(r"$\mathit \eta_{{\mathrm{{bf}}}}$")

        axFFvel.set_xlim(0,1)
        axFFvel.set_ylim(0,1)

        plt.savefig(doFF.paperFigurePath+"{}compareBFFs.pdf".format(doFF.codeRawLocalCO12 ), bbox_inches='tight'  )


    def reLabelData(self,dataLabel,TB):

        """
        remove cloud not in TB, used to remove bad channels,finish this today
        :param data:
        :param TB:
        :return:
        """

        #coValues=dataLabel[dataLabel>0]
        labelSets= doFF.getLabelSet(dataLabel) #[Z0, Y0, X0, clusterValue1D ]

        newLabel=np.zeros_like(dataLabel)

        for eachCloud in TB:
            cloudID= eachCloud["_idx"]
            cloudIndices=doFF.getIndices(labelSets,cloudID)

            newLabel[cloudIndices] = cloudID

        return newLabel





    def getCuoffListAndError(self,calCode, rawCOFITS,labelCOFITS,rmsFITS,saveTag="", reCalculate=True ):

        """

        :param rawCOFITS:
        :param labelCOFITS:
        :return:
        """
        xSave=  self.saveOverallCutoff+saveTag+"xCutoff"
        ySave=  self.saveOverallCutoff+saveTag+"yCutoff"
        yErrSave=  self.saveOverallCutoff+saveTag+"yErrCutoff"

        #need to remove bad channels




        if not reCalculate and  os.path.isfile(xSave+".npy") and  os.path.isfile(ySave+".npy") and  os.path.isfile(yErrSave+".npy"):
            x = np.load(xSave + ".npy" )
            yList = np.load(ySave + ".npy" )
            yErrList = np.load(yErrSave + ".npy" )

            #print x
            #print yList
            #print yErrList , "???????????????????????"
            #return x,yList, None

            return x,yList, yErrList


        FFTB=self.getFFTB(calCode)

        TB=Table.read( FFTB )
        dataCO,headCO=myFITS.readFITS(rawCOFITS)
        labelCO,_ = myFITS.readFITS( labelCOFITS )

        rmsData,rmsHead = myFITS.readFITS(rmsFITS)

        x =  doFF.cutoffFactors

        #need to remove bad channels


        labelCO=self.reLabelData( labelCO,TB  )



        coValues=dataCO[labelCO>0]
        tmbInRMS = dataCO/rmsData

        rmsValue=tmbInRMS[labelCO>0]

 
        yList=[]

        yErrList = []



        for eachX  in x:
            subCO= coValues[rmsValue>=eachX ]

            y=np.sum(subCO)*doFF.getVelResolution()
            yErr= np.sqrt(  len(subCO) )*doFF.getMeanRMS()* doFF.getVelResolution()
            if yErr==0:
                yErr=   doFF.getMeanRMS()* doFF.getVelResolution()

            yErrList.append(  yErr )

            yList.append(y)




        np.save(xSave,x)
        np.save(ySave,yList)
        np.save(yErrSave,yErrList)


        
        return x,yList, yErrList




    def getOverallCutoff(self,calCode,calID=0,reCalculate=False,dim=5):
        """

        :param calCode:
        :return:
        """

        #getRaw
        doFF.calCode=calCode

        rawCO=doFF.getRawCOFITS()
        labelFITS= doFF.getSmoothAndNoiseFITSSingle( smFactor=1.0,  noiseFactor=0.0 , getCleanFITS=True)
        rmsFITS=doFF.getRMSFITS()


        #need to remove bad sources



        x,y,yError= self.getCuoffListAndError(calCode, rawCO,labelFITS,rmsFITS,saveTag=calCode,reCalculate=reCalculate)


        # drawTheFigure
        fig = plt.figure(figsize=(10, 8))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 22, 'serif': ['Helvetica']})

        if 1:
            mpl.rcParams['text.latex.preamble'] = [
                r'\usepackage{tgheros}',  # helvetica font
                r'\usepackage{sansmath}',  # math-font matching  helvetica
                r'\sansmath',  # actually tell tex to use it!
                r'\usepackage{siunitx}',  # micro symbols
                r'\sisetup{detect-all}',  # force siunitx to use the fonts
            ]

        axFitting = fig.add_subplot(1, 1, 1)
        # axFitting.scatter(x, y, s=15, color='red', label="The noise RMS is 0.5 K")
        # axFitting.scatter(x, x*y, s=15, color='red', label="The noise RMS is 0.5 K")
        # axFitting.scatter(  x  ,  y , s=15,  color='red'   ) #
        #axFitting.errorbar(x, -y, yerr=yError, markersize=2, fmt='o', c='red', capsize=0.5, elinewidth=1, lw=1)

        axFitting.errorbar(  x,   y ,   markersize=2, fmt='o', c='red', capsize=0.5, elinewidth=1, lw=1)








        ##fitting with polynomial


        fitingX=x

        fittingy=y
        fittingYerror = yError



        #print doFF.fittingCutOFF(fitingX, fittingy, fittingYerror, dim=  dim )
        params, paramas_covariance, useFunction ,fitData= doFF.fittingCutOFF(fitingX, fittingy, fittingYerror, dim= dim )
        label = ""
        if 1:
            mwispFF, mwispFFError = doFF.getFFandErrorCutoff( params, paramas_covariance,  targetSNR=2 ,dimension= len(params)  )
            if np.isnan( mwispFF ) or np.isinf(mwispFF):
                mwispFF=0
                mwispFFError=0


            label= "$\eta_{{ \\rm sen}}={:.3f}\pm{:.3f}$".format(mwispFF ,mwispFFError  ) ##"Quadratic with Gaussian CDF "

        if 0: #test

            polyX = x[x<=mu]  # x[0:8]
            polyY =  y[x<=mu]  # y[0:8]
            axDiff.scatter(polyX, polyY, color='blue')

            p0=np.polyfit(polyX,polyY,deg= 2 )

            p = np.poly1d(p0)

            axDiff.plot(polyX,  p(polyX),  c='black' , lw=1)
            #mwispFF, mwispFFError = self.getFFandErrorCutoff( params, paramas_covariance,  targetSNR=2 ,dimension= len(params)  )

        fittingX = np.arange(0, np.max(fitData[0] ), 0.01)
        # axFitting.plot( fittingX  ,  useFunction(fittingX,params[0], params[1] , params[2]   ), color='blue'  )
        # if useCutoff:
        # fittingX = np.arange( - np.max(x), np.max(x), 0.01)

        #useFunction= ffFunctionCutTest


        #if useFunction == ffFunctionCutTest:
            #if np.isnan(mwispFFError) or np.isinf(mwispFFError):
                #mwispFFError=999


        if 0:
             # with astropy
            fitter = modeling.fitting.LevMarLSQFitter()
            model1 = modeling.models.Gaussian1D(amplitude=np.mean(y), mean=0, stddev=5)
            model1.mean.fixed = True

            model2 = modeling.models.Gaussian1D(amplitude=np.mean(y) / 2, mean=0, stddev=3)
            model2.mean.fixed = True

            gInit = model1   + model2

            fitted_model = fitter(gInit, x, y, weights=1. / yError)
            # fitted_model = fitter(gInit, x, y)
            # print gInit.mean.value , gInit.stddev.value ,gInit.amplitude.value

            # print fitted_model.mean_0.value,  fitted_model.stddev_0.value, fitted_model.amplitude_0.value
            # print fitted_model.mean_1.value,  fitted_model.stddev_1.value, fitted_model.amplitude_1.value
            # print model2.mean.value ,model2.stddev.value ,model2.amplitude.value

            ffs = axFitting.plot(fittingX, fitted_model(fittingX), color='blue', lw=1.0  )
        else:
            #p  ass
            pass
            ffs = axFitting.plot(fittingX, useFunction(fittingX, *params), color='blue', lw=1.0 ,label=label     )


        axFitting.axvline(x=0, ls="--", color='black',lw=1.2,alpha=0.5)
        armStr = doFF.getArmStr()
        lineStr = doFF.getLineStr()
        fillingWMISP = r"{} molecular clouds in the {} arm".format(lineStr, armStr)
        ffs2 = axFitting.plot(fittingX, useFunction(fittingX, *params), color='white', alpha=0, lw=1.0,
                              label=fillingWMISP)

        axFitting.legend(loc=1, handlelength=1, fontsize=20)

        axFitting.set_xlabel(r"Cutoff ($\sigma$)")
        # axFitting.set_ylabel("Flux (arcmin)")
        axFitting.set_ylabel(r"Flux ($\rm K\ km\ s$$^{-1}$ $\mathrm{\Omega}_\mathrm{A}$)")
        ###########################
        saveTagFigue = "{}_factorFitting_ID{}{}_dim{}".format(doFF.calCode, calID, "",dim)
        #axFitting.set_xlim(0,8)

        if calID != 0:

            if individuals and not testPoly:
                saveTagFigue = saveTagFigue + "_indivi"
                plt.savefig(self.figurePath + "{}.pdf".format(saveTagFigue), bbox_inches='tight')


            if testPoly:
                saveTagFigue = saveTagFigue + "_indivi"
                plt.savefig(self.figurePath + "{}Poly.pdf".format(saveTagFigue), bbox_inches='tight')

            if not testPoly and not individuals:
                plt.savefig(self.figurePath + "{}.png".format(saveTagFigue), bbox_inches='tight', dpi=100)



        else:
            plt.savefig(doFF.paperFigurePath + "{}.png".format(saveTagFigue), bbox_inches='tight', dpi=100)
            plt.savefig(doFF.paperFigurePath + "{}.pdf".format(saveTagFigue), bbox_inches='tight')

        plt.close(fig)
        gc.collect()

        return  mwispFF, mwispFFError

        print  useFunction(2 , *params) /  useFunction(0, *params) , "BFFFF of ",calCode
        #return mwispFF, 0, fittingParaAndError, paramas_covariance

    def zzz(self):
        """

        :return:
        """
        pass



#\eta_{\mathrm{bf}}

doMain=fillingMain()


###################
if 1:
    #doFF.calCode=doFF.codeRawLocalCO12 #codeRawLocalCO12
    doFF.calCode=doFF.codeRawLocalCO12 #codeRawLocalCO12

    #doFF.callFillingFactorAllSM( ) #this is over all cloud fits
    doFF.printFillingCat()
    #doFF.printFluxCat()


    sys.exit()

if 0:
    doMain.drawCutoffBFF(doFF.codeRawLocalCO13 , dim=  5 ,showX2=True   )


    doMain.drawCutoffBFF(doFF.codeRawLocalCO12 , dim=  5 ,showX2=True   )

    #doMain.drawCutoffBFF(doFF.codeRawLocalCO12 , dim=  5    )



if 0:


    for eachCode in doFF.allRawCodeList:

        doMain.getCutffFF(  eachCode   , drawFigure=False  , cutByPeak=False ,dim= 5 )
        #doMain.getCutffFF(  eachCode   , drawFigure=False , cutByPeak=False ,dim=3 )



    sys.exit()


    doFF.calCode = doFF.codeRawLocalCO12
    TB = Table.read(  "cutoffFlux_pureEdgeInfo_fillingFactor_fluxTB_rawLocalCO12rawLocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit")

    TB.sort("sum")
    TB.reverse()
    doFF.calculateFillingFactorCutoff(TB, drawFigure=False ,  dim=5)
    #doFF.calculateFillingFactorCutoff(TB, drawFigure=True, inputID=213561, saveTag="cutoffBFFcase", dim=5)
    #doFF.calculateFillingFactorCutoff(TB, drawFigure=True, inputID=113544, saveTag="cutoffBFFcase", dim=5)








if 0:

    for eachCode in doFF.allRawCodeList:
        doMain.getOverallCutoff( eachCode  ,reCalculate=  False   ,dim= 5     ) #can get error

    #doMain.getOverallCutoff( doFF.codeRawLocalCO12 ,reCalculate=  False     )
    #doMain.getOverallCutoff( doFF.codeRawLocalCO12  ,reCalculate=  False   ,dim= 3    )
    #doMain.getOverallCutoff( doFF.codeRawLocalCO12  ,reCalculate=  False   ,dim= 5    )

    #doMain.getOverallCutoff(doFF.codeRawSgrCO12 ,reCalculate=  False     )
    #doMain.getOverallCutoff(doFF.codeRawSgrCO13 ,reCalculate=  False     )

    #doMain.getOverallCutoff(doFF.codeRawLocalCO13  ,reCalculate=  True       )
    #doMain.getOverallCutoff(doFF.codeRawLocalCO18  ,reCalculate=  False     )
    sys.exit()

if 0:
    # doMain.getCuoffFlux( doFF.codeRawLocalCO18 )
    # doMain.getCutffFF( doFF.codeRawLocalCO13 )

    doFF.calCode = doFF.codeRawLocalCO13
    TB = Table.read( "cutoffFlux_5pureEdgeInfo_fillingFactor_fluxTB_rawLocalCO13rawLocalCO13_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit")

    TB.sort("sum")
    TB.reverse()
    #doFF.calculateFillingFactorCutoff(TB, drawFigure=True, inputID=228690, saveTag="cutoffBFFcase", dim=5)
    doFF.calculateFillingFactorCutoff(TB, drawFigure=True, inputID=24687, saveTag="cutoffBFFcase", dim=5)
    doFF.calculateFillingFactorCutoff(TB, drawFigure=True, inputID=37758, saveTag="cutoffBFFcase", dim=5)

    doFF.calculateFillingFactorCutoff(TB, drawFigure=True,dim= 5 )
    sys.exit()

if 0:
    # doMain.getCuoffFlux( doFF.codeRawLocalCO18 )
    # doMain.getCutffFF( doFF.codeRawLocalCO13 )

    doFF.calCode = doFF.codeRawLocalCO12
    TB = Table.read(
        "cutoffFlux3_pureEdgeInfo_fillingFactor_fluxTB_rawLocalCO12rawLocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit")

    TB.sort("sum")
    TB.reverse()

    #doFF.calculateFillingFactorCutoff(TB, drawFigure=True, inputID=434975, saveTag="cutoffBFFcase", dim=5)
    #doFF.calculateFillingFactorCutoff(TB, drawFigure=True, inputID=392595, saveTag="cutoffBFFcase", dim=5)
    #doFF.calculateFillingFactorCutoff(TB, drawFigure=True, inputID=  430167 , saveTag="cutoffBFFcase", dim=5)
    #doFF.calculateFillingFactorCutoff(TB, drawFigure=True, inputID=  260664 , saveTag="cutoffBFFcase", dim=5)
    doFF.calculateFillingFactorCutoff(TB, drawFigure=True, inputID=  462967 , saveTag="cutoffBFFcase", dim= 5 )

    doFF.calculateFillingFactorCutoff(TB, drawFigure=True, inputID=  332289 , saveTag="cutoffBFFcase", dim=5 )

    #doFF.calculateFillingFactorCutoff(TB, drawFigure=True)
    sys.exit()




if 0:
     #doFF.drawCloudSizeChange(drawCode="min")
    #doFF.drawCloudSizeChange(drawCode="max")
    #doFF.drawCloudSizeChange(drawCode="mean")
    doFF.calCode=doFF.codeRawLocalCO12
    doFF.drawCloudNumberChange()

    doFF.calCode=doFF.codeRawLocalCO13
    doFF.drawCloudNumberChange()

    doFF.calCode=doFF.codeRawLocalCO18
    doFF.drawCloudNumberChange()





if 0: #draw figures for sensitivities
    pass
    #figure 1, two cases of cutoff specific molecular clouds
    #figure 3, relationship between the two BFFs

    #figure4, overall fiting of all clouds
    #display the results of surveys, in total five figures

    #pput those figures in discussion

    #first, display two cases

    #draw first Table
    doFF.calCode=doFF.codeRawLocalCO12
    fluxTB=Table.read( "rawLocalCO12BFFclean.fit" )
    if 1:
        doFF.calculateFillingFactorCutoff( fluxTB,drawFigure=True,inputID=369770,saveTag="cutoffBFFcase", dim=5)
        doFF.calculateFillingFactorCutoff( fluxTB,drawFigure=True,inputID=344936,saveTag="cutoffBFFcase", dim=5)
        doFF.calculateFillingFactorCutoff( fluxTB,drawFigure=True,inputID=323137,saveTag="cutoffBFFcase", dim=5)

    #relation between the two BFFs
    if 1:
        doMain.drawBFFrelatin( fluxTB )

    if 1:
        combinTB=doMain.getCombinCutoffTable()

        doMain.drawCutoffBFF( doFF.codeRawLocalCO12, saveTag="cutoffBFFAllMWISP",inputTB=combinTB ,errorThresh=0.2, showLegend=True,cloudSource="MWISP molecular clouds")


    if 1:
        doMain.drawCutOFFBFFSurvey()

    sys.exit()



if 0:

    #for eachCode in doFF.allRawCodeList:
    for eachCode in doFF.surveyCodeList:
        doMain.drawCutoffBFF( eachCode    )





if 0:
    doFF.calCode= doFF.codeRawLocalCO12
    #doMain.produceSameSizeClouds()
    doMain.compareCloudSameSize()





if  0:
    doMain.fittingFFWithAllTBs(doFF.allRawCodeList,"AllClouds")
    #doMain.fittingFFWithAllTBs(doFF.rawOutCodeList,"AllOutClouds")
    #doMain.fittingFFWithAllTBs(doFF.rawLocalCodeList,"AllLocalClouds")
    #doMain.fittingFFWithAllTBs(doFF.rawSgrCodeList,"smgrClouds")
    #doMain.fittingFFWithAllTBs(doFF.rawScuCodeList,"AllScuClouds")

    CO12All=[doFF.codeRawLocalCO12,doFF.codeRawSgrCO12,doFF.codeRawScuCO12,doFF.codeRawOutCO12]
    CO13All=[doFF.codeRawLocalCO13,doFF.codeRawSgrCO13,doFF.codeRawScuCO13,doFF.codeRawOutCO13]
    CO18All=[doFF.codeRawLocalCO18,doFF.codeRawSgrCO18,doFF.codeRawScuCO18,doFF.codeRawOutCO18]

    #doMain.fittingFFWithAllTBs( CO12All ,"AllCO12Clouds")
    #doMain.fittingFFWithAllTBs( CO13All ,"AllCO13Clouds")
    #doMain.fittingFFWithAllTBs( CO18All ,"AllCO18Clouds")

    sys.exit()



if  0: # draw survey results

    doFF.drawSurveyBFF()
    sys.exit()



if 0: #testing

    #for eachCode in doFF.allRawCodeList:
    for eachCode in [ doFF.codeRawLocalCO12]:


        doFF.calCode=  eachCode  #doFF.codeRawLocalCO18

        #doFF.getFluxListForEachCloud(calCode=  eachCode  )

        doMain.getFillingFactorAndEdge(eachCode, drawFigure=True)

    sys.exit()



if 0: #single test of draw filing factor


    for i in [1,0]:
        if i: #draw individual
            testPoly = False
            drawIndi = True
        else: #draw polynomial
            testPoly=True
            drawIndi=True

        doMain.drawFillingFactor(doFF.codeRawLocalCO12,drawFigure=True,inputID=95454, testPoly=testPoly, individuals=drawIndi )
        doMain.drawFillingFactor(doFF.codeRawLocalCO12,drawFigure=True,inputID=95336, testPoly=testPoly, individuals=drawIndi )
        doMain.drawFillingFactor(doFF.codeRawLocalCO12,drawFigure=True,inputID=115887, testPoly=testPoly, individuals=drawIndi )
        doMain.drawFillingFactor(doFF.codeRawLocalCO12,drawFigure=True,inputID=129123, testPoly=testPoly, individuals=drawIndi )

    sys.exit()





if 0:
    doMain.drawChiSquareTest(trainingRatio=0.8)
    doMain.drawChiSquareTest(trainingRatio=0.7)

    sys.exit()




if 0: #pipeline of calculating filling factors


    allCodeList =  doFF.allRawCodeList+doFF.surveyCodeList

    for processCode in  allCodeList:

        if processCode in   doFF.rawLocalCodeList  or processCode in  doFF.rawSgrCodeList :
            continue



        doFF.calCode=processCode
        # SMOOTH
        if 0:
            doFF.calCode= processCode
            doFF.smoothFITSbySMFactor()
        # addnoise

        if 1:
            smFiles = doFF.getSmFITSFileList()

            for eachSMF in smFiles:
                print "Processing ", eachSMF
                doFF.addNoiseByRMSFITS(eachSMF, noiseFactor=0.0)
        #search cloud
        if 1:
            doMain.cleanFITS(processCode, onlyFirst=False)
            doMain.removeUselessFITS(processCode)

        #calculate BFF
        if 1:
            doFF.getFluxListForEachCloud(calCode=  processCode  )
        if 1:
            doMain.getFillingFactorAndEdge(  processCode, drawFigure=False )



    sys.exit()









if 0:

    doMain.testPeakBFF() #examin










if 0:
    doFF.drawBFFDiscuss()
    sys.exit()



if 0:
    doFF.calCode=doFF.codeRawLocalCO12 #codeRawLocalCO12

    #doFF.printCloudNumberTable()
    doFF.printFillingCat()
    #doFF.printFluxCat()


    sys.exit()






if  0: #testing

    for eachCode in doFF.allRawCodeList:
    #for eachCode in [ doFF.codeRawLocalCO12,doFF.codeRawSgrCO12,doFF.codeRawScuCO12,doFF.codeRawOutCO12   ]:
    #for eachCode in [ doFF.codeRawLocalCO13,doFF.codeRawSgrCO13,doFF.codeRawScuCO13,doFF.codeRawOutCO13   ]:
    #for eachCode in [ doFF.codeRawLocalCO18,doFF.codeRawSgrCO18,doFF.codeRawScuCO18,doFF.codeRawOutCO18   ]:

    #for eachCode in doFF.allRawCodeList :

        doFF.calCode=  eachCode  #doFF.codeRawLocalCO18

        tb=doMain.getFFTB(doFF.calCode)
        testTB=Table.read(tb)

        doFF.fittingAngularSizeFF( tb,showSizeRange=[-1,50],useODR=False ,showColor="velocity")

    sys.exit()








if 0: #testing

    for eachCode in doFF.allRawCodeList:
    #for eachCode in [ doFF.codeRawLocalCO12,doFF.codeRawSgrCO12,doFF.codeRawScuCO12,doFF.codeRawOutCO12   ]:
    #for eachCode in [ doFF.codeRawLocalCO13,doFF.codeRawSgrCO13,doFF.codeRawScuCO13,doFF.codeRawOutCO13   ]:
    #for eachCode in [ doFF.codeRawLocalCO18,doFF.codeRawSgrCO18,doFF.codeRawScuCO18,doFF.codeRawOutCO18   ]:

    #for eachCode in doFF.allRawCodeList :

        doFF.calCode=  eachCode  #doFF.codeRawLocalCO18

        tb=doMain.getFFTB(doFF.calCode)
        testTB=Table.read(tb)

        doFF.fittingAngularSizeFF( tb,showSizeRange=[-1,50],useODR=False )

    sys.exit()





if 0: #pipeline of calculating filling factors

    processCode=   doFF.codeRawLocalCO13
    doFF.calCode=processCode

    if 1:
        doMain.cleanFITS(processCode, onlyFirst=False)
        #doMain.removeUselessFITS(processCode)


    sys.exit()

if 0:

    doFF.calCode=doFF.codeRawLocalCO13
    doFF.testEqualCutOff()
    #doFF.testCloudVoxelRatio()



    sys.exit()


if 0:
    doFF.calCode=doFF.surveyCodeCfACO12
    doFF.getSNRMedian()

    doFF.calCode=doFF.surveyCodeGRSCO13
    doFF.getSNRMedian()

    doFF.calCode=doFF.surveyCodeOGSCO12
    doFF.getSNRMedian()

    doFF.calCode=doFF.surveyCodeCOHRSCO32
    doFF.getSNRMedian()
    sys.exit()






if 0: #prepare the data
    #doFF.calCode=  doFF.surveyCodeCfACO12

    doFF.clearSurveyData()





if 0: #

    processCode=   doFF.surveyCodeOGSCO12
    # SMOOTH
    doFF.calCode= processCode

    doMain.drawFillingFactor( processCode ,drawFigure=True,inputID=  26  )

    sys.exit()



if 0: #pipeline

    processCode=   doFF.surveyCodeCfACO12
    # SMOOTH
    doFF.calCode= processCode


    #calculate BFF
    doFF.getFluxListForEachCloud(calCode=  processCode  )
    doMain.getFillingFactorAndEdge(  processCode, drawFigure=False )

    sys.exit()














if 0: #pipeline

    processCode=   doFF.surveyCodeOGSCO12
    # SMOOTH
    doFF.calCode= processCode
    doFF.smoothFITSbySMFactor()
    # addnoise

    if 1:
        smFiles = doFF.getSmFITSFileList()

        for eachSMF in smFiles:
            print "Processing ", eachSMF
            doFF.addNoiseByRMSFITS(eachSMF, noiseFactor=0.0)
    #search cloud
    doMain.cleanFITS(processCode, onlyFirst=False)
    doMain.removeUselessFITS(processCode)

    #calculate BFF
    doFF.getFluxListForEachCloud(calCode=  processCode  )
    doMain.getFillingFactorAndEdge(  processCode, drawFigure=False )



    sys.exit()






if 0: #pipeline

    processCode=   doFF.codeRawLocalCO12
    # SMOOTH
    doFF.calCode= processCode

    doMain.drawFillingFactor( processCode ,drawFigure=True,inputID= 485654  )

    sys.exit()

if 0: #pipeline

    processCode=   doFF.surveyCodeGRSCO13
    # SMOOTH
    doFF.calCode= processCode


    #calculate BFF
    #doFF.getFluxListForEachCloud(calCode=  processCode  )
    #doMain.getFillingFactorAndEdge(  processCode, drawFigure=True,inputID=)
    doMain.drawFillingFactor( processCode ,drawFigure=True,inputID= 26779  )



    sys.exit()













if 0:

    for eachCode in doFF.allRawCodeList :

        doMain.getFillingFactorAndEdge(  eachCode, drawFigure=False )

    sys.exit()




if 0:  # addNoise
    # codeList=[doFF.codeRawSgrCO13,doFF.codeRawSgrCO18,doFF.codeRawScuCO12,doFF.codeRawScuCO13, doFF.codeRawScuCO18,doFF.codeRawOutCO12,doFF.codeRawOutCO13,doFF.codeRawOutCO18  ]
    # for eachCode in  codeList:
    eachCode = doFF.codeRawLocalCO13
    doFF.calCode = eachCode

    if 1:
        smFiles = doFF.getSmFITSFileList()

        for eachSMF in smFiles:
            print "Processing ", eachSMF
            doFF.addNoiseByRMSFITS(eachSMF, noiseFactor=0.0)

    doMain.cleanFITS(eachCode, onlyFirst=False)
    doMain.removeUselessFITS(eachCode)

    sys.exit()




if 0:  # find clouds #

    for eachCode in doFF.allRawCodeList:

        if eachCode==doFF.codeRawLocalCO12:
            continue

        doFF.calCode=eachCode
        doFF.getFluxListForEachCloud(calCode=  eachCode  )

    sys.exit()













if 0: # addNoise
    codeList=[doFF.codeRawSgrCO13,doFF.codeRawSgrCO18,doFF.codeRawScuCO12,doFF.codeRawScuCO13, doFF.codeRawScuCO18,doFF.codeRawOutCO12,doFF.codeRawOutCO13,doFF.codeRawOutCO18  ]
    for eachCode in  codeList:
        eachCode= eachCode
        doFF.calCode= eachCode

        if 1:
            smFiles = doFF.getSmFITSFileList()

            for eachSMF in smFiles:
                print "Processing ",eachSMF
                doFF.addNoiseByRMSFITS( eachSMF,noiseFactor=0.0 )

        doMain.cleanFITS( eachCode  , onlyFirst=False )
        doMain.removeUselessFITS( eachCode  )
    sys.exit()









if 0: # clean all fits, #Do this tonight

    #for eachCode in doFF.allRawCodeList: #dow SgrCO12 and Local CO12 first

    eachCode=doFF.codeRawSgrCO12
    doMain.cleanFITS( eachCode  , onlyFirst=False )
    doMain.removeUselessFITS( eachCode  )

    sys.exit()







if 0: # clean all fits, #Do this tonight

    #for eachCode in doFF.allRawCodeList: #dow SgrCO12 and Local CO12 first

    eachCode=doFF.codeRawLocalCO12
    doMain.cleanFITS( eachCode  , onlyFirst=False )
    doMain.removeUselessFITS( eachCode  )

    sys.exit()



if 0: # addNoise
    doFF.calCode=doFF.codeRawLocalCO12
    smFiles = doFF.getSmFITSFileList()

    for eachSMF in smFiles:
        print "Processing ",eachSMF
        doFF.addNoiseByRMSFITS( eachSMF,noiseFactor=0.0 )






if 0:

    for eachCode in doFF.allRawCodeList:
        doFF.calCode=  eachCode
        doFF.smoothFITSbySMFactor()

    sys.exit()



if 0: #clean data
    #some area has too large noise, need to removethis par
    pass
    doFF.cleanRawData()









#the first flux,should used the raw molecular flux, not 2.6 sigma cut, fix this but






if 0: # clean all fits, #Do this tonight

    for eachCode in doFF.allRawCodeList:

        doMain.cleanFITS( eachCode, onlyFirst=False )

        doMain.removeUselessFITS(eachCode)

    sys.exit()






if 0:

    for eachCode in doFF.allRawCodeList :

        doFF.calCode= eachCode  #doFF.codeRawLocalCO12
        tbFile=doMain.getFFTB(doFF.codeRawLocalCO12)

        doFF.calculateFillingFactor( tbFile, drawFigure=True )
        #doFF.calculateFillingFactor( tbFile,   inputID=395765 )
        aaaa
    sys.exit()


if 0:
    doFF.calCode=doFF.codeRawLocalCO12

    TB=Table.read("cutoffFF_pureEdgeInfo_fillingFactor_fluxTB_rawLocalCO12rawLocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit")


    newTB=doFF.addMWISPFFerrorCutoff(TB)
    newTB= doFITS.selectTBByColRange(newTB, "cutFFfillingFactorMWISPError",maxV=0.3  )


    newTB.write("aaaaa.fit")

    sys.exit()














if 0: #testing

    #for eachCoe in doFF.allRawCodeList:
    for eachCode in [ doFF.codeRawLocalCO12]:
    #for eachCode in doFF.allRawCodeList :

        doFF.calCode=  eachCode  #doFF.codeRawLocalCO18

        tb=doMain.getFFTB(doFF.calCode)
        testTB=Table.read(tb)

        doFF.testThreeFunctions( tb  )

    sys.exit()











if 0:#find clouds #
    doFF.calCode=doFF.codeRawLocalCO18
    noiseFiles=doFF.getSmoothAndNoiseFITSSingle( smFactor=1.0, noiseFactor=0.0 )

    print noiseFiles
    doFF.cleanFITSsigma2( noiseFiles  )

    sys.exit()









if 0:

    doMain.getFillingFactorAndEdge(doFF.codeRawScuCO12,drawFigure=False)

    sys.exit()





#draw filling factor for a particular ID




if 0:

    for eachCode in doFF.allRawCodeList :
        doMain.getFillingFactorAndEdge(eachCode,drawFigure=False)



    sys.exit()


if 0: # clean all fits, #Do this tonight

    eachCode=doFF.codeRawSgrCO12
    doMain.cleanFITS( eachCode, onlyFirst=False )
    doMain.removeUselessFITS(eachCode)

    sys.exit()


if 0:  # find clouds #

    #for eachCode in doFF.allRawCodeList:
    doFF.getFluxListForEachCloud(calCode= doFF.codeRawSgrCO12 )

    sys.exit()

#########

if 0: # addNoise
    doFF.calCode=doFF.codeRawSgrCO12
    smFiles = doFF.getSmFITSFileList()

    for eachSMF in smFiles:
        print "Processing ",eachSMF
        doFF.addNoiseByRMSFITS( eachSMF,noiseFactor=0.0 )



if 0: # test remove bad TB
    doFF.calCode=doFF.codeRawOutCO18

    #TBTest="fluxTB_rawScuCO12rawScuCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    #TBTest="fluxTB_rawSgrCO18rawSgrCO18_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    #TBTest="fluxTB_rawScuCO12rawScuCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    TBTest="fluxTB_rawOutCO18rawOutCO18_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"

    tb=Table.read( TBTest )

    cloudLeft=doFF.removeFakeClouds(tb) #removeByLBVrange(tb,vRange=[0,30],bRange=[-5,5],lRange=[30,40])
    cloudLeft.write("cloudLeftTest.fit",overwrite=True)

    sys.exit()



if 0:
    doFF.removeBadSgrChannelCO12()

if 0:
    for i in [26,27, 28,29,30,31 ]:
        doMain.testBackChannelFFT(i)



if 0:  # find clouds #


    for eachCode in doFF.allRawCodeList:
        doFF.getFluxListForEachCloud(calCode= eachCode )

    sys.exit()

if 0: # test remove bad TB
    doFF.calCode=doFF.codeRawOutCO12

    TBTest="fillingFactor_fluxTB_rawOutCO12rawOutCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"

    tb=Table.read( TBTest )

    cloudLeft=doFF.removeFakeClouds(tb) #removeByLBVrange(tb,vRange=[0,30],bRange=[-5,5],lRange=[30,40])
    cloudLeft.write("cloudLeftTest.fit",overwrite=True)

    sys.exit()



if 0: # clean all fits, #Do this tonight

    for eachCode in doFF.allRawCodeList:

        doMain.cleanFITS( eachCode, onlyFirst=False )
        doMain.removeUselessFITS(eachCode)

    sys.exit()



#calculating flux Scu
if 0:
    doFF.calCode = doFF.codeRawOutCO13
    doFF.getFluxListForEachCloud()



    sys.exit()





if 0:#find clouds #
    doFF.calCode=doFF.codeRawOutCO13
    noiseFiles=doFF.getSmoothListFixNoise(noiseFactor=0.0 )

    for eachF in noiseFiles:
        doFF.cleanFITSsigma2( eachF )

    sys.exit()

if 0: #thest how to remove fake clouds
    #doMain.getFillingFactorAndEdge( doFF.codeRawOutCO12 )
    #doMain.getFillingFactorAndEdge( doFF.codeRawLocalCO12 )
    #doMain.getFillingFactorAndEdge( doFF.codeRawSgrCO12)

    #doMain.getFillingFactorAndEdge( doFF.codeRawLocalCO12)
    #doMain.getFillingFactorAndEdge( doFF.codeRawOutCO18 )
    doMain.getFillingFactorAndEdge( doFF.codeRawOutCO13 )

    sys.exit()


if 0:
    doFF.calCode=doFF.codeRawScuCO18

    lRange,bRange=doFITS.box(43.5674609,-3.7843809,1157.026 ,680.262 ,0)

    showChannel = 144
    doFF.testBadChannel( showChannel,lRange,bRange)
    sys.exit()



if 0:

    #for eachCode in doFF.allRawCodeList:
    doMain.getFillingFactorAndEdge(eachCode)

    sys.exit()




if 0:
    for eachCode in doFF.allRawCodeList:
        doFF.calCode= eachCode
        doFF.produceMaskSigma2(cutoffSigma=2)
    sys.exit()








if 0:
    doFF.calCode=doFF.codeRawLocalCO12
    doFF.addAllNoise()
    sys.exit()


if 0:

    doFF.calCode=doFF.codeRawLocalCO12
    doFF.getFluxListForEachCloudNoiseChange()
    sys.exit()







if 0:  #drawLocalCO13


    #codeList1=[doFF.codeRawLocalCO12, doFF.codeRawLocalCO13, doFF.codeRawLocalCO18 ]
    codeList1=[doFF.codeRawOutCO12, doFF.codeRawOutCO13, doFF.codeRawOutCO18 ]

    #codeList2=[doFF.codeRawScuCO12,doFF.codeRawScuCO13,doFF.codeRawScuCO18]
    codeList2=[ ]

    codeListAll=codeList1+codeList2

    for eachCode in codeListAll:
        doFF.calCode = eachCode
        drawTB = doMain.getFFTB( doFF.calCode )

        #doFF.drawFillingRelation(  drawTB, drawCode=doFF.drawCodeFlux)
        #doFF.drawFillingRelation(  drawTB, drawCode=doFF.drawCodeArea)
        doFF.drawFillingRelation( drawTB, drawCode=doFF.drawCodeSize)

    sys.exit()






if 0:

    for eachCode in doFF.allRawCodeList:
        doFF.calCode= eachCode
        tb=doMain.getFFTB(doFF.calCode)
        doFF.fittingAngularSizeFF( tb )
    sys.exit()

if 0: #check particular clouds, smaller filling factors
    doFF.calCode=doFF.codeRawLocalCO12

    checkTB="edgeInfo_fillingFactor_fluxTB_rawLocalCO12rawLocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    cleanTB=Table.read(checkTB)


    #getTB= doFF.selectBySizeAndFactorRange(checkTB,sizeRange=[0, 4 ],factorRange=[0.7, 1])

    getTB= doFF.selectBySizeAndFactorRange(checkTB,sizeRange=[5, 10 ],factorRange=[0.0, 0.1])

    getTB=doFF.addMWISPFFerror(getTB)

    getTB.sort(doFF.ffMWISPErrorCol)
    #print getTB[doFF.ffMWISPCol]

    printTB= doFITS.selectTBByCols(getTB,["_idx","area_exact","sum","v_cen","peak","pixN",doFF.ffMWISPCol,doFF.ffMWISPErrorCol])
    printTB.sort(doFF.ffMWISPErrorCol)

    print printTB
    print len(printTB)
    #testIDList=[2607,10959,6656,117419,29631]
    #testIDList=[274907, 215799]
    #aaaaaaaa
    for eachRow in getTB:

        #doFF.calFFByID(testCode,eachID ,cleanTB )
        eachID=eachRow["_idx"]
        doFF.drawInMap(doFF.calCode,eachID)
        doFF.calFFByRow(eachRow,drawFigure=True)

    sys.exit()


    #draw those figures







if 0:  #drawLocalCO13


    codeList1=[doFF.codeRawLocalCO12, doFF.codeRawLocalCO13, doFF.codeRawLocalCO18 ]
    #codeList2=[doFF.codeRawScuCO12,doFF.codeRawScuCO13,doFF.codeRawScuCO18]
    codeList2=[ ]

    codeListAll=codeList1+codeList2

    for eachCode in codeListAll:
        doFF.calCode = eachCode
        drawTB = doMain.getFFTB( doFF.calCode )

        #doFF.drawFillingRelation(  drawTB, drawCode=doFF.drawCodeFlux)
        #doFF.drawFillingRelation(  drawTB, drawCode=doFF.drawCodeArea)
        doFF.drawFillingRelation( drawTB, drawCode=doFF.drawCodeSize)

    sys.exit()








if 0: #clean the first fits, to do the flux calculation

    for eachCode in doFF.allRawCodeList:

        doMain.cleanFITS( eachCode, onlyFirst=True )



if 0:
    doMain.removeUselessFITS(doFF.codeRawScuCO18)



if 0:
    doFF.calCode=doFF.codeRawLocalCO12 #usually only the first cloud would be examined
    doFF.cloudStat()
    doFF.calCode=doFF.codeRawLocalCO13 #usually only the first cloud would be examined
    doFF.cloudStat()





if 0: #test distance layers
    doFF.calCode=doFF.codeRawSgrCO12
    doFF.testDistanceLayers()
    #doFF.produceLVByRemovingTheLargestCloud()
    sys.exit()


if 0:  #drawLocalCO13

    for eachCode in [doFF.codeRawSgrCO12]:
        doFF.calCode = eachCode
        drawTB = doMain.getFFTB( doFF.calCode )

        #doFF.drawFillingRelation(  drawTB, drawCode=doFF.drawCodeFlux)
        #doFF.drawFillingRelation(  drawTB, drawCode=doFF.drawCodeArea)
        doFF.drawFillingRelation( drawTB, drawCode=doFF.drawCodeSize)



if 0:

    doMain.getFillingFactorAndEdge(doFF.codeRawLocalCO12)
    doMain.getFillingFactorAndEdge(doFF.codeRawLocalCO13)
    doMain.getFillingFactorAndEdge(doFF.codeRawLocalCO18)

    doMain.getFillingFactorAndEdge(doFF.codeRawOutCO12)
    doMain.getFillingFactorAndEdge(doFF.codeRawOutCO13)
    doMain.getFillingFactorAndEdge(doFF.codeRawOutCO18)


    sys.exit()
    #

if 0:

    doMain.getFillingFactorAndEdge(doFF.codeRawScuCO12)
    doMain.getFillingFactorAndEdge(doFF.codeRawScuCO13)
    doMain.getFillingFactorAndEdge(doFF.codeRawScuCO18)

    sys.exit()

if 0:

    doMain.getFillingFactorAndEdge(doFF.codeRawSgrCO12)
    doMain.getFillingFactorAndEdge(doFF.codeRawSgrCO13)
    doMain.getFillingFactorAndEdge(doFF.codeRawSgrCO18)

    sys.exit()

#calculating flux Scu
if 0:
    doFF.calCode = doFF.codeRawScuCO12
    doFF.getFluxListForEachCloud()

    doFF.calCode = doFF.codeRawScuCO13
    doFF.getFluxListForEachCloud()

    doFF.calCode = doFF.codeRawScuCO18
    doFF.getFluxListForEachCloud()



    sys.exit()





if 0: ### #used to test which clouds are at the edge of data cube
    doFF.calCode=doFF.codeRawLocalCO12
    cleanFITSName= doFF.getSmoothAndNoiseFITSSingle(getCleanFITS=True)  #doFF.tmpPath+  #"LocalCO12LocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fits"
    print cleanFITSName
    tbName =  "fillingFactor_fluxTB_rawLocalCO12rawLocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    doFF.getCloudCubes( cleanFITSName ,  tbName,writeFITS=True)
    sys.exit()


if 0: #get Sgr

    doFF.calCode = doFF.codeRawSgrCO12
    doFF.getFluxListForEachCloud()

    doFF.calCode = doFF.codeRawSgrCO13
    doFF.getFluxListForEachCloud()

    doFF.calCode = doFF.codeRawSgrCO18
    doFF.getFluxListForEachCloud()


    sys.exit()





if 0:


    file2="edgeInfo_fluxTB_rawOutCO12rawOutCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    doFF.calCode=doFF.codeRawOutCO12

    doFF.calculateFillingFactor(file2, drawFigure=True)


    sys.exit()






if 0:  #drawLocalCO13

    drawTB = "edgeInfo_fillingFactor_fluxTB_rawLocalCO12rawLocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    doFF.calCode = doFF.codeRawLocalCO12

    #doFF.drawFillingRelation(  drawTB, drawCode=doFF.drawCodeFlux)

    #doFF.drawFillingRelation(  drawTB, drawCode=doFF.drawCodeArea)
    doFF.drawFillingRelation( drawTB, drawCode=doFF.drawCodeSize)



    drawTB = "fillingFactor_edgeInfo_fluxTB_rawLocalCO13rawLocalCO13_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    doFF.calCode = doFF.codeRawLocalCO13

    #doFF.drawFillingRelation(  drawTB, drawCode=doFF.drawCodeFlux)

    #doFF.drawFillingRelation(  drawTB, drawCode=doFF.drawCodeArea)
    doFF.drawFillingRelation( drawTB, drawCode=doFF.drawCodeSize)


    drawTB = "fillingFactor_edgeInfo_fluxTB_rawLocalCO18rawLocalCO18_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    doFF.calCode = doFF.codeRawLocalCO18

    #doFF.drawFillingRelation(doFF.codeLocalCO12, drawTB, drawCode=doFF.drawCodeFlux)

    #doFF.drawFillingRelation(doFF.codeLocalCO12, drawTB, drawCode=doFF.drawCodeArea)
    doFF.drawFillingRelation( drawTB, drawCode=doFF.drawCodeSize)


    drawTB = "fillingFactor_edgeInfo_fluxTB_rawOutCO12rawOutCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    doFF.calCode = doFF.codeRawOutCO12

    #doFF.drawFillingRelation(doFF.codeLocalCO12, drawTB, drawCode=doFF.drawCodeFlux)

    #doFF.drawFillingRelation(doFF.codeLocalCO12, drawTB, drawCode=doFF.drawCodeArea)
    doFF.drawFillingRelation( drawTB, drawCode=doFF.drawCodeSize)



    drawTB = "fillingFactor_edgeInfo_fluxTB_rawOutCO13rawOutCO13_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    doFF.calCode = doFF.codeRawOutCO13

    #doFF.drawFillingRelation(doFF.codeLocalCO12, drawTB, drawCode=doFF.drawCodeFlux)

    #doFF.drawFillingRelation(doFF.codeLocalCO12, drawTB, drawCode=doFF.drawCodeArea)
    doFF.drawFillingRelation( drawTB, drawCode=doFF.drawCodeSize)

    drawTB = "fillingFactor_edgeInfo_fluxTB_rawOutCO18rawOutCO18_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    doFF.calCode = doFF.codeRawOutCO18

    #doFF.drawFillingRelation(doFF.codeLocalCO12, drawTB, drawCode=doFF.drawCodeFlux)

    #doFF.drawFillingRelation(doFF.codeLocalCO12, drawTB, drawCode=doFF.drawCodeArea)
    doFF.drawFillingRelation( drawTB, drawCode=doFF.drawCodeSize)





    sys.exit()



if 0:




    file1="edgeInfo_fluxTB_rawLocalCO18rawLocalCO18_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    file2="edgeInfo_fluxTB_rawOutCO12rawOutCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    file3="edgeInfo_fluxTB_rawOutCO13rawOutCO13_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    file4="edgeInfo_fluxTB_rawOutCO18rawOutCO18_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    file5="edgeInfo_fluxTB_rawLocalCO13rawLocalCO13_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"




    doFF.calculateFillingFactor(file1, drawFigure=False)
    doFF.calculateFillingFactor(file2, drawFigure=False)
    doFF.calculateFillingFactor(file3, drawFigure=False)
    doFF.calculateFillingFactor(file4, drawFigure=False)
    doFF.calculateFillingFactor(file5, drawFigure=False)

    sys.exit()




if 0: ### #used to test which clouds are at the edge of data cube
    doFF.calCode=doFF.codeRawLocalCO13
    cleanFITSName= doFF.getSmoothAndNoiseFITSSingle(getCleanFITS=True)  #doFF.tmpPath+  #"LocalCO12LocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fits"
    print cleanFITSName
    tbName =  "fluxTB_rawLocalCO13rawLocalCO13_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    doFF.getCloudCubes( cleanFITSName ,  tbName,writeFITS=False)


    doFF.calCode=doFF.codeRawLocalCO18
    cleanFITSName= doFF.getSmoothAndNoiseFITSSingle(getCleanFITS=True)  #doFF.tmpPath+  #"LocalCO12LocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fits"
    print cleanFITSName
    tbName =  "fluxTB_rawLocalCO18rawLocalCO18_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    doFF.getCloudCubes(  cleanFITSName ,  tbName,writeFITS=False)




    doFF.calCode=doFF.codeRawOutCO12
    cleanFITSName= doFF.getSmoothAndNoiseFITSSingle(getCleanFITS=True)  #doFF.tmpPath+  #"LocalCO12LocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fits"
    print cleanFITSName
    tbName =  "fluxTB_rawOutCO12rawOutCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    doFF.getCloudCubes(  cleanFITSName ,  tbName,writeFITS=False)


    doFF.calCode=doFF.codeRawOutCO13
    cleanFITSName= doFF.getSmoothAndNoiseFITSSingle(getCleanFITS=True)  #doFF.tmpPath+  #"LocalCO12LocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fits"
    print cleanFITSName
    tbName =  "fluxTB_rawOutCO13rawOutCO13_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    doFF.getCloudCubes(  cleanFITSName ,  tbName,writeFITS=False)


    doFF.calCode=doFF.codeRawOutCO18
    cleanFITSName= doFF.getSmoothAndNoiseFITSSingle(getCleanFITS=True)  #doFF.tmpPath+  #"LocalCO12LocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fits"
    print cleanFITSName
    tbName =  "fluxTB_rawOutCO18rawOutCO18_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    doFF.getCloudCubes(  cleanFITSName ,  tbName,writeFITS=False)






if 0: #get flux tb

    doFF.calCode = doFF.codeRawLocalCO13
    doFF.getFluxListForEachCloud()

    doFF.calCode = doFF.codeRawLocalCO18
    doFF.getFluxListForEachCloud()

    doFF.calCode = doFF.codeRawOutCO12
    doFF.getFluxListForEachCloud()

    doFF.calCode = doFF.codeRawOutCO13
    doFF.getFluxListForEachCloud()

    doFF.calCode = doFF.codeRawOutCO18
    doFF.getFluxListForEachCloud()

    sys.exit()


#########

if 0: # addNoise
    doFF.calCode=doFF.codeRawOutCO18
    smFiles = doFF.getSmFITSFileList()

    for eachSMF in smFiles:
        print "Processing ",eachSMF
        doFF.addNoiseByRMSFITS( eachSMF,noiseFactor=0.0 )


if 0:#find clouds #
    doFF.calCode=doFF.codeRawOutCO13
    noiseFiles=doFF.getSmoothListFixNoise(noiseFactor=0.0 )
    for eachF in noiseFiles:
        doFF.cleanFITSsigma2( eachF )


if 0:#find clouds #
    doFF.calCode=doFF.codeRawOutCO18
    noiseFiles=doFF.getSmoothListFixNoise(noiseFactor=0.0 )
    for eachF in noiseFiles:
        doFF.cleanFITSsigma2( eachF )


#####################

if 0:#find clouds #
    doFF.calCode=doFF.codeRawLocalCO18
    noiseFiles=doFF.getSmoothListFixNoise(noiseFactor=0.0 )
    for eachF in noiseFiles:
        doFF.cleanFITSsigma2( eachF )
    sys.exit()


if 0:  # raw velocity resolutuion

    # doFF.calCode= doFF.codeOutCO12
    # doFF.smoothFITSbySMFactor(doFF.outCO12FITS)
    # sys.exit()

    ##### Local
    #doFF.calCode= doFF.codeRawLocalCO12
    #doFF.smoothFITSbySMFactor(doFF.rawLocalCO12FITS)

    #doFF.calCode= doFF.codeRawLocalCO13
    #doFF.smoothFITSbySMFactor(doFF.rawLocalCO13FITS)

    #doFF.calCode= doFF.codeRawLocalCO18
    #doFF.smoothFITSbySMFactor(doFF.rawLocalCO18FITS)

    #out
    #doFF.calCode= doFF.codeRawOutCO12
    #doFF.smoothFITSbySMFactor(doFF.rawOutCO12FITS)

    #doFF.calCode = doFF.codeRawOutCO13
    #doFF.smoothFITSbySMFactor(doFF.rawOutCO13FITS)

    #doFF.calCode = doFF.codeRawOutCO18
    #doFF.smoothFITSbySMFactor(doFF.rawOutCO18FITS)


    #raw Sgr
    doFF.calCode= doFF.codeRawSgrCO12
    doFF.smoothFITSbySMFactor(doFF.rawSgrCO12FITS)

    doFF.calCode= doFF.codeRawSgrCO13
    doFF.smoothFITSbySMFactor(doFF.rawSgrCO13FITS)

    doFF.calCode= doFF.codeRawSgrCO18
    doFF.smoothFITSbySMFactor(doFF.rawSgrCO18FITS)

    #raw Scu
    doFF.calCode= doFF.codeRawScuCO12
    doFF.smoothFITSbySMFactor(doFF.rawScuCO12FITS)

    doFF.calCode= doFF.codeRawScuCO13
    doFF.smoothFITSbySMFactor(doFF.rawScuCO13FITS)

    doFF.calCode= doFF.codeRawScuCO18
    doFF.smoothFITSbySMFactor(doFF.rawScuCO18FITS)

    sys.exit()






if 0:  # calculate filling factors for each molecular clouds

    drawTB = "edgeInfo_fillingFactor_fluxTB_rawLocalCO12rawLocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    doFF.calCode = doFF.codeRawLocalCO12

    #doFF.drawFillingRelation(doFF.codeLocalCO12, drawTB, drawCode=doFF.drawCodeFlux)

    #doFF.drawFillingRelation(doFF.codeLocalCO12, drawTB, drawCode=doFF.drawCodeArea)
    doFF.drawFillingRelation(doFF.codeLocalCO12, drawTB, drawCode=doFF.drawCodeSize)

    sys.exit()

if 0: ### #used to test which clouds are at the edge of data cube
    doFF.calCode=doFF.codeRawLocalCO12
    cleanFITSName= doFF.getSmoothAndNoiseFITSSingle(getCleanFITS=True)  #doFF.tmpPath+  #"LocalCO12LocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fits"
    print cleanFITSName
    tbName =  "fillingFactor_fluxTB_rawLocalCO12rawLocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    doFF.getCloudCubes(doFF.codeRawLocalCO12,doFF.rawLocalCO12FITS, cleanFITSName ,  tbName,writeFITS=False)



    #calTBFile="fluxTB_rawLocalCO12rawLocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"
    #doFF.calculateFillingFactor( calTBFile ,drawFigure=False)

    sys.exit()




if 0:
    doFF.calCode=doFF.codeRawLocalCO12

    doFF.getFluxListForEachCloud()

    sys.exit()





if 0:
    doFF.checkRMS()
    #doFF.checkData()


if 0:  # out CO 12

    # doFF.calCode= doFF.codeOutCO12
    # doFF.smoothFITSbySMFactor(doFF.outCO12FITS)
    # sys.exit()

    ##### Local
    #doFF.calCode= doFF.codeRawLocalCO12
    #doFF.smoothFITSbySMFactor(doFF.rawLocalCO12FITS)

    doFF.calCode= doFF.codeRawLocalCO13
    doFF.smoothFITSbySMFactor(doFF.rawLocalCO13FITS)

    doFF.calCode= doFF.codeRawLocalCO18
    doFF.smoothFITSbySMFactor(doFF.rawLocalCO18FITS)

    #out
    #doFF.calCode= doFF.codeRawOutCO12
    #doFF.smoothFITSbySMFactor(doFF.rawOutCO12FITS)

    doFF.calCode = doFF.codeRawOutCO13
    doFF.smoothFITSbySMFactor(doFF.rawOutCO13FITS)

    doFF.calCode = doFF.codeRawOutCO18
    doFF.smoothFITSbySMFactor(doFF.rawOutCO18FITS)


    sys.exit()











if 0:#find clouds #
    doFF.calCode=doFF.codeLocalCO18
    noiseFiles=doFF.getSmoothListFixNoise(noiseFactor=0.0 )
    for eachF in noiseFiles:
        doFF.cleanFITSsigma2( eachF )
    sys.exit()
if 0:
    doFF.calCode=doFF.codeLocalCO12
    fileList=doFF.getSmoothListFixNoise(noiseFactor=1.0)

    print fileList


if 0:  # out CO 12

    # doFF.calCode= doFF.codeOutCO12
    # doFF.smoothFITSbySMFactor(doFF.outCO12FITS)
    # sys.exit()

    ##### Local

    doFF.calCode = doFF.codeLocalCO12
    doFF.smoothFITSbySMFactor(doFF.localCO12FITS)

    doFF.calCode = doFF.codeLocalCO13
    doFF.smoothFITSbySMFactor(doFF.localCO13FITS)

    doFF.calCode = doFF.codeLocalCO18
    doFF.smoothFITSbySMFactor(doFF.localCO18FITS)

    ################
    doFF.calCode = doFF.codeOutCO12
    doFF.smoothFITSbySMFactor(doFF.outCO12FITS)


    doFF.calCode = doFF.codeOutCO13
    doFF.smoothFITSbySMFactor(doFF.outCO13FITS)

    doFF.calCode = doFF.codeOutCO18
    doFF.smoothFITSbySMFactor(doFF.outCO18FITS)

    ##############Sgr
    doFF.calCode = doFF.codeSgrCO12
    doFF.smoothFITSbySMFactor(doFF.sgrCO12FITS)

    doFF.calCode = doFF.codeSgrCO13
    doFF.smoothFITSbySMFactor(doFF.sgrCO13FITS)

    doFF.calCode = doFF.codeSgrCO18
    doFF.smoothFITSbySMFactor(doFF.sgrCO18FITS)

    sys.exit()

    ##### Scu
    doFF.calCode = doFF.codeScuCO12
    doFF.smoothFITSbySMFactor(doFF.scuCO12FITS)

    doFF.calCode = doFF.codeScuCO13
    doFF.smoothFITSbySMFactor(doFF.scuCO13FITS)

    doFF.calCode = doFF.codeScuCO18
    doFF.smoothFITSbySMFactor(doFF.scuCO18FITS)







if 0:
    doFF.produceAllRMSFITS()


if 0: #check particular clouds, smaller filling factors

    checkTB="LocalCO12FillingFactorTBAll.fit"

    testCode=doFF.codeLocalCO12

    getTB= doFF.selectBySizeAndFactorRange(checkTB,sizeRange=[0, 4 ],factorRange=[0.75,1])
    getTB.sort(doFF.ffMWISPCol)

    print getTB["area_exact"]
    print getTB["sum"]
    print getTB["pixN"]
    print getTB["peak"]



    testIDList=[2607,10959,6656,117419,29631]
    for eachID in testIDList:

        doFF.calFFByID(testCode,eachID,drawFigure=True,useSigmaCut=True)
        doFF.drawInMap(testCode,eachID)


    sys.exit()


if 0: #check particular clouds, smaller filling factors

    checkTB="LocalCO12FillingFactorTBAll.fit"

    testCode=doFF.codeLocalCO12
    testID= 27090

    #getTB= doFF.selectBySizeAndFactorRange(checkTB,sizeRange=[5,7.5],factorRange=[0,0.3])
    #getTB.sort("area_exact")

    #print getTB

    doFF.calFFByID(testCode,testID,drawFigure=True,useSigmaCut=True)
    doFF.drawInMap(testCode,testID)


    sys.exit()



#get systematic parameters and run the script

if  0: #calculate filling factors for each molecular clouds
    #to get enough memory

    calCode= sys.argv[1]


    #if calCode==doFF.codeLocalCO12:
        #sys.exit()

    doFF.getFFForEachCloud( calCode, drawFigure=True, useSigmaCut=True, calAllCloud=True)
    sys.exit()


if 0:
    tbName =  "LocalCO12FillingFactorTBAll.fit"
    doFF.getFluxListForEachCloud( tbName )






if 0:
    #pass
    #recalculate fits and keep the dbscan and clean fits only for local fits

    for eachCode, eachN in zip(doFF.scuCodeList, doFF.noiseList  ):
        doMain.cleanByCode(eachCode,eachN,removeFITS=False)


    for eachCode, eachN in zip(doFF.outCodeList, doFF.noiseList  ):
        doMain.cleanByCode(eachCode,eachN,removeFITS=False)



    #sys.exit()




if 0:
    #pass
    #recalculate fits and keep the dbscan and clean fits only for local fits

    for eachCode, eachN in zip(doFF.localCodeList, doFF.noiseList  ):
        doMain.cleanByCode(eachCode,eachN,removeFITS=False)

    for eachCode, eachN in zip(doFF.sgrCodeList, doFF.noiseList  ):
        doMain.cleanByCode(eachCode,eachN,removeFITS=False)

    #sys.exit()

if 0:
    #doFF.printFillingCat()

    sys.exit()

if 0: #calculate and save

    doFF.recordFillingFactorData()
    sys.exit()



if 0: #test CO18

    #doMain.drawFillingFactor(doFF.codeLocalCO12,doFF.MWISPrmsCO12)
    #doMain.drawFillingFactor(doFF.codeLocalCO13,doFF.MWISPrmsCO13)
    #doMain.drawFillingFactor(doFF.codeLocalCO18,doFF.MWISPrmsCO18)

    #doMain.drawFillingFactor(doFF.codeOutCO13, doFF.MWISPrmsCO13)

    doMain.drawFillingFactor(doFF.codeOutCO13, doFF.MWISPrmsCO13)

    sys.exit()


if 0: #test Main
    #doMain.addNoiseToSmFiles(doFF.codeOutCO12,doFF.MWISPrmsCO12)

    #doMain.cleanByCode(doFF.codeOutCO12,doFF.MWISPrmsCO12)
    #doMain.drawFillingFactor(doFF.codeOutCO12,doFF.MWISPrmsCO12)


    #doMain.addNoiseToSmFiles(doFF.codeOutCO13,doFF.MWISPrmsCO13)
    #doMain.cleanByCode(doFF.codeOutCO13, doFF.MWISPrmsCO13)
    #doMain.drawFillingFactor(doFF.codeOutCO13,doFF.MWISPrmsCO13)

    codeList=[doFF.codeLocalCO12, doFF.codeOutCO12, doFF.codeSgrCO12, doFF.codeScuCO12, \
              doFF.codeLocalCO13, doFF.codeOutCO13, doFF.codeSgrCO13, doFF.codeScuCO13, \
              doFF.codeLocalCO18, doFF.codeOutCO18, doFF.codeSgrCO18, doFF.codeScuCO18,]



    noiseList=[doFF.MWISPrmsCO12, doFF.MWISPrmsCO12, doFF.MWISPrmsCO12, doFF.MWISPrmsCO12, \
              doFF.MWISPrmsCO13, doFF.MWISPrmsCO13,  doFF.MWISPrmsCO13, doFF.MWISPrmsCO13, \
              doFF.MWISPrmsCO18, doFF.MWISPrmsCO18,  doFF.MWISPrmsCO18, doFF.MWISPrmsCO18,]


    for i in range(len(codeList)):
        code=codeList[i]
        print "Calculating code ", code

        noise=noiseList[i]

        doMain.drawFillingFactor(code,noise)
        #doMain.addNoiseToSmFiles( code ,noise )
        #doMain.cleanByCode(code, noise )
    sys.exit()

if 0:
    doFF.getArmSubFTS() #produce






if 0:


    doFF.calCode="Q1LocalCO13" #-6-30 km/s
    #doFF.calCode="Q1Sgr" #30-70 km/s

    doFF.drawFluxOfNoise()




if 0: #clean Q1LocalCO13

    doFF.calCode="Q1LocalCO13" #30-70 km/s

    smFiles= doFF.getSmFITSFileList()

    for eachF in smFiles:
        doFF.addNoiseSingle(eachF,absoluteNoise=doFF.MWISPrmsCO12 )

    Q1LocalCO13SmFITS = doFF.getAbsNoiseFITSName()

    for eachF in Q1LocalCO13SmFITS:
        #continue  # Do not do this a second time

        doFF.cleanFITSsigma2(eachF, removeFITS=False)



if 0:#clean Sgr

    doFF.calCode="Q1Sgr" #30-70 km/s

    Q1SgrSmFITS = doFF.getAbsNoiseFITSName()

    print Q1SgrSmFITS

    for eachF in Q1SgrSmFITS:
        #continue  # Do not do this a second time

        doFF.cleanFITSsigma2(eachF, removeFITS=False)

if 0:# local 13CO

    doFF.calCode="Q1LocalCO13" #30-70 km/s
    doFF.rawFITS="./data/G2650LocalCO13.fits"

    #smFiles= doFF.getSmFITSFileList()

    #add noise to 0.5 K

    #for eachF in smFiles:
        #doFF.addNoiseSingle(eachF,absoluteNoise=doFF.MWISPrmsCO12 )


    #step 1 smooth by factor
    doFF.smoothFITSbySMFactor(doFF.rawFITS)

    #step 2, add noise to all smoothe filts
    #doFF.addNoseToAllSmFiles( doClean=True)

    sys.exit()





if 0:
    doFF.calCode = "Q1Local"  # deal with local
    doFF.drawRMSwithBEAM()



if 0:#pipeline line

    doFF.calCode="Q1Sgr" #30-70 km/s
    doFF.rawFITS="./data/G2650Sgr70.fits"

    smFiles= doFF.getSmFITSFileList()

    #add noise to 0.5 K

    for eachF in smFiles:
        doFF.addNoiseSingle(eachF,absoluteNoise=doFF.MWISPrmsCO12 )


    #step 1 smooth by factor
    #doFF.smoothFITSbySMFactor(doFF.rawFITS)

    #step 2, add noise to all smoothe filts
    #doFF.addNoseToAllSmFiles( doClean=True)



if 0:
    doFF.calCode = "Q1Local"  # deal with local
    doFF.factorFitting()


if 0:
    doFF.calCode = "Q1Local"  # deal with local

    doFF.pipeLineTestabsK()


if 0:#test sensitivity

    doFF.calCode = "Q1Local"
    testFile=   "/home/qzyan/WORK/myDownloads/fillingFactor/tmpFiles/Q1LocalG2650Local30_SmFactor_2.0_noise_0.5absK.fits"

    doFF.compareSensitivity(testFile)

if 0:
    #doFF.calCode="Q1Local"

    #testFile=   "/home/qzyan/WORK/myDownloads/fillingFactor/tmpFiles/Q1LocalG2650Local30_SmFactor_2.0_noise_0.5absK.fits"
    #testTBFile= "/home/qzyan/WORK/myDownloads/fillingFactor/tmpFiles/Q1LocalG2650Local30_SmFactor_2.0_noise_0.5absKdbscanS2P4Con1_Clean.fit"

    doFF.getEqualCutOff( testFile,  testTBFile  )




if 0:
    doFF.calCode="Q1Local"

    doFF.addNoiseToWMISPNoise()

if 0: #draw table
    doFF.calCode="Q1Local"

    #tbList,NoiseArray=doFF.getTBFixSmFactor()
    #print doFF.getTotalFlux(tbList)
    doFF.drawBeamEffect()
    doFF.drawNoiseEffect()


if 0:# test DBSCAN
    #cleanFITSsigma2
    doFF.calCode="smallTest"
    doFF.rawFITS="./data/subTest.fits" #the raw CO fits is used to calculate the total flux
    doFF.cleanFITSsigma2(os.path.join(doFF.tmpPath,"smallTestsubTest_SmFactor_4.0.fits"), removeFITS=True)



if 0:#pipeline Q1Local
    pass
    doFF.calCode="Q1Local"
    doFF.rawFITS="./data/G2650Local30.fits.fits"

    #step 1 smooth by factor
    #doFF.smoothFITSbySMFactor(doFF.rawFITS)

    #step 2, add noise to all smoothe filts
    doFF.addNoseToAllSmFiles( doClean=True)



if 0:  # pipeline test

    doFF.calCode="smallTest"
    doFF.rawFITS="./data/subTest.fits"

    # step 1 smooth by factor
    #doFF.smoothFITSbySMFactor(doFF.rawFITS)

    # step 1, add noise to all smoothe filts
    doFF.addNoseToAllSmFiles(doClean=True)


if 0:#pipeline test

    doFF.calCode="pipeLineTest"
    doFF.rawFITS="./data/subTest.fits"

    #step 1 smooth by factor
    doFF.smoothFITSbySMFactor(doFF.rawFITS)

    #step 1, add noise to all smoothe filts
    doFF.addNoseToAllSmFiles()




if 0:

    fitsList= doFF.getFITSnamesByPrefix(doFF.tmpPath,"subTest_SM_")
    doFF.calRMS(fitsList)

    #doFF.calRMS("./data/G2650Local30.fits")


if 0:#test pipline

    rawCOFITS="./data/subTest.fits"

    doFF.smoothFITSbySMFactor(rawCOFITS)
if 0:
    doFF.drawfluxChange()


if 0:
    doFF.smoothCube( doFF.CO12Raw ,disFactor=2 )
    doFF.smoothCube( doFF.CO12Raw ,disFactor=3 )
    doFF.smoothCube( doFF.CO12Raw ,disFactor=4 )
    doFF.smoothCube( doFF.CO12Raw ,disFactor=5 )
    doFF.smoothCube( doFF.CO12Raw ,disFactor=6 )
    doFF.smoothCube( doFF.CO12Raw ,disFactor=7 )
    doFF.smoothCube( doFF.CO12Raw ,disFactor=8 )

if 0:


    for i in [2,3,4,5,6,7,8]:
        smFITS= "smoothedCO12Factor{}.fits".format(i)
        doFF.addNoise( smFITS,0.5 )





#doFF.addNoise( "resampledCO12.fits",0.5/4, 0.5  )