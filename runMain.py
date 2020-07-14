#from fillingfactor import checkFillingFactor
from fillingfactor import *

import os
import sys
from astropy.table import Table,vstack
from myPYTHON import *

doFITS=myFITS() #used to deal fits with myPYTHON

doFF=checkFillingFactor()



class fillingMain(object):

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
            from matplotlib.colors import LogNorm
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

    def fittingFFWithAllTBs(self,codeList,saveTag=""):
        """
        fitting an overall filling factor, with all samples
        :return:
        """

        #according to the code list

        #for eachCoe in doFF.allRawCodeList:

        combineTB=None

        for eachCode in  codeList:
        #for eachCode in doFF.allRawCodeList :

            doFF.calCode=  eachCode  #doFF.codeRawLocalCO18
            tb=doMain.getFFTB(doFF.calCode)
            testTB=Table.read(tb)
            print "Total number of clouds of ",eachCode, " is ",len(testTB)
            #add size in beam colName



            testTB[doFF.sizeInBeamCol] =doFF.getCloudSize(testTB)/doFF.getBeamSize()

            if combineTB is None:
                combineTB=testTB

            else:
                combineTB=vstack([combineTB,testTB])



        # use BeamFactorUnitFor combined clouds
        doFF.fittingAngularSizeFF( combineTB,showSizeRange=[-1,40],saveTag=saveTag,useBeamFactorUnit=True )



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

        tb = self.getFFTB(doFF.calCode)
        ffTB = Table.read(tb)
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


        for eachEcut in  errorCutList:

            processTB=doFF.selectByErrorCut(ffTB,errorCut=eachEcut)

            trainingTB, testTB = doFF.splitTBIntoTrainingTest(processTB, trainingRatio=trainingRatio)



            if 1:

                cs1,para1,paraError1=doFF.testThreeFunctionsByTraining(testffAndSizeFunc1,   trainingTB,testTB, trainingRatio=trainingRatio, errorCut = eachEcut  )
                chiSquareList1.append(cs1)

                cs2,para2,paraError2=doFF.testThreeFunctionsByTraining(testffAndSizeFunc2,   trainingTB,testTB,  trainingRatio=trainingRatio,  errorCut = eachEcut  )
                chiSquareList2.append(cs2)


                cs3,para3,paraError3 =doFF.testThreeFunctionsByTraining(testffAndSizeFunc3,  trainingTB,testTB,  trainingRatio=trainingRatio,  errorCut = eachEcut  )
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
        rc('font', **{'family': 'sans-serif', 'size': 13, 'serif': ['Helvetica']})
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]


        axTest = plt.subplot(122 )

        axTest.plot(  errorCutList, chiSquareList1,'.-',color='blue', markersize=7,lw=1,label=r"$\mathit f=\mathit l/\left(\mathit l+ a\right)$")
        axTest.plot(  errorCutList, chiSquareList2,'.-',color='green', markersize=7,lw=1,label=r"$\mathit f=a\left(1-\exp\left(-b\mathit l\right)\right)$")
        axTest.plot(  errorCutList, chiSquareList3,'.-',color='red', markersize=7,lw=1,label=r"$\mathit f= a\exp\left(-b \mathit l\right)  +c$")

        axTest.legend(loc=4)
        axTest.set_ylim(0,0.15)

        axTest.set_xlabel("Maximum relative error of filling factors")
        axTest.set_ylabel("Chi-Square")
        axTest.set_ylabel("Weighted rms of test data")

        #print "{:.0f}\% test data".format(1-trainingRatio), "What the hell?????????"
        at = AnchoredText(r"{:.0f}\% test data".format((1-trainingRatio)*100), loc=2, frameon=False)
        axTest.add_artist(at)

        saveTag= "modelTestTestRatio{:.2f}".format(1-trainingRatio)


        #draw a figure with
        axFF = plt.subplot(121 )

        print "Traing number of examples",len(trainHalfErrorCut)
        print "test number   of examples",len(testHalfErrorCut)


        drawXTrain = doFF.getCloudSize(trainHalfErrorCut)
        drawYTrain = trainHalfErrorCut[doFF.ffMWISPCol]
        yErrorTrain = trainHalfErrorCut[doFF.ffMWISPErrorCol]

        axFF.errorbar(drawXTrain,   drawYTrain ,  yerr= yErrorTrain ,    markersize=1 , fmt='o',  color='gray' , capsize=0.1,  elinewidth=0.5 , lw=1 ,alpha=0.8  ,label="Training data (50\% maximum relative BFF errors)" )


        drawXTest = doFF.getCloudSize(testHalfErrorCut)
        drawYTest = testHalfErrorCut[doFF.ffMWISPCol]
        yErrorTest = testHalfErrorCut[doFF.ffMWISPErrorCol]

        axFF.errorbar(drawXTest,   drawYTest ,  yerr= yErrorTest ,    markersize=1 , fmt='o',  color='black' , capsize=0.1,  elinewidth=0.5 , lw=1 ,alpha=0.8  ,label="Validation data (50\% maximum relative BFF errors)" )
        showSizeRange=[-1,40]
        axFF.set_ylim([-0.01, 1.1])
        axFF.set_xlim( showSizeRange )

        axFF.legend(loc= 4 )

        #draw four lines
        x = np.arange(-0.1, showSizeRange[1], 0.01)
        axFF.plot(x, testffAndSizeFunc1(x, *paraHalfErrorCutff1 ), "b-", lw=1 )
        axFF.plot(x, testffAndSizeFunc2(x, *paraHalfErrorCutff2 ), "g-", lw=1 )
        axFF.plot(x, testffAndSizeFunc3(x, *paraHalfErrorCutff3 ), "r-", lw=1 )

        axFF.axvline(x=0, ls="--", color='black', lw=0.8)
        axFF.axhline(y=0,ls="--",color='black',lw=0.8)

        plt.savefig(   saveTag+".png"  , bbox_inches='tight', dpi=600)
        plt.savefig(   saveTag+".pdf"  , bbox_inches='tight' )







    def getCutffFF(self,calCode):

        doFF.calCode=calCode

        rawCOFITS=doFF.getRawCOFITS()
        labelFITS= doFF.getSmoothAndNoiseFITSSingle( smFactor=1.0,  noiseFactor=0.0 , getCleanFITS=True)
        TBFile=self.getFFTB(calCode)
        rmsFITS=doFF.getRMSFITS()
        doFF.getSmoothFluxColCutoff(rawCOFITS,labelFITS,rmsFITS, TBFile )




    def zzz(self):
        """

        :return:
        """
        pass





doMain=fillingMain()

if 0: #prepare the data
    #doFF.calCode=  doFF.surveyCodeCfACO12

    doFF.clearSurveyData()


if 0: #testing

    for eachCode in doFF.allRawCodeList:
    #for eachCode in [ doFF.codeRawLocalCO12]:
    #for eachCode in doFF.allRawCodeList :

        doFF.calCode=  eachCode  #doFF.codeRawLocalCO18

        tb=doMain.getFFTB(doFF.calCode)
        testTB=Table.read(tb)

        doFF.fittingAngularSizeFF( tb,showSizeRange=[-1,40],useODR=False )

    sys.exit()




if 0:
    doMain.fittingFFWithAllTBs(doFF.allRawCodeList,"AllClouds")
    #doMain.fittingFFWithAllTBs(doFF.rawOutCodeList,"AllOutClouds")
    #doMain.fittingFFWithAllTBs(doFF.rawLocalCodeList,"AllLocalClouds")
    #doMain.fittingFFWithAllTBs(doFF.rawSgrCodeList,"AllSgrClouds")
    #doMain.fittingFFWithAllTBs(doFF.rawScuCodeList,"AllScuClouds")

    CO12All=[doFF.codeRawLocalCO12,doFF.codeRawSgrCO12,doFF.codeRawScuCO12,doFF.codeRawOutCO12]
    CO13All=[doFF.codeRawLocalCO13,doFF.codeRawSgrCO13,doFF.codeRawScuCO13,doFF.codeRawOutCO13]
    CO18All=[doFF.codeRawLocalCO18,doFF.codeRawSgrCO18,doFF.codeRawScuCO18,doFF.codeRawOutCO18]

    #doMain.fittingFFWithAllTBs( CO12All ,"AllCO12Clouds")
    #doMain.fittingFFWithAllTBs( CO13All ,"AllCO13Clouds")
    #doMain.fittingFFWithAllTBs( CO18All ,"AllCO18Clouds")

    sys.exit()




if 0: #

    processCode=   doFF.surveyCodeOGSCO12
    # SMOOTH
    doFF.calCode= processCode

    doMain.drawFillingFactor( processCode ,drawFigure=True,inputID=  26  )

    sys.exit()



if 0: #pipeline

    processCode=   doFF.surveyCodeMWISPCO13
    # SMOOTH
    if 1:
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

if 0: # draw survey results
    processCode=   doFF.surveyCodeCfACO12
    doFF.calCode=  doFF.surveyCodeCfACO12

    doFF.drawSurveyBFF()
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







if   0:
    doFF.calCode=doFF.codeRawLocalCO12 #codeRawLocalCO12

    #for eachCode in [doFF.codeRawLocalCO12, doFF.codeRawLocalCO13, doFF.codeRawLocalCO18]:
    #for eachCode in [doFF.codeRawOutCO12 ]:

        #doFF.calCode=eachCode
        #doFF.drawCloudNumberChange()

    doFF.callFillingFactorAllSM( ) #this is over all cloud fits
    #doFF.printFillingCat()
    #doFF.printFluxCat()


    sys.exit()


if 0:
    doFF.drawBFFDiscuss()
    sys.exit()

if 0: #testing

    #for eachCoe in doFF.allRawCodeList:
    #for eachCode in [ doFF.codeRawLocalCO12,doFF.codeRawSgrCO12,doFF.codeRawScuCO12,doFF.codeRawOutCO12   ]:
    #for eachCode in [ doFF.codeRawLocalCO13,doFF.codeRawSgrCO13,doFF.codeRawScuCO13,doFF.codeRawOutCO13   ]:
    for eachCode in [ doFF.codeRawLocalCO18,doFF.codeRawSgrCO18,doFF.codeRawScuCO18,doFF.codeRawOutCO18   ]:

    #for eachCode in doFF.allRawCodeList :

        doFF.calCode=  eachCode  #doFF.codeRawLocalCO18

        tb=doMain.getFFTB(doFF.calCode)
        testTB=Table.read(tb)

        doFF.fittingAngularSizeFF( tb,showSizeRange=[-1,40],useODR=False )

    sys.exit()





if 0: #single test of draw filing factor
    if 1: #draw individual
        testPoly = False
        drawIndi = True
    else: #draw polynomial
        testPoly=True
        drawIndi=True

    doMain.drawFillingFactor(doFF.codeRawLocalCO12,drawFigure=True,inputID=466713, testPoly=testPoly, individuals=drawIndi )
    doMain.drawFillingFactor(doFF.codeRawLocalCO12,drawFigure=True,inputID=488352, testPoly=testPoly, individuals=drawIndi )
    doMain.drawFillingFactor(doFF.codeRawLocalCO12,drawFigure=True,inputID=465842, testPoly=testPoly, individuals=drawIndi )
    doMain.drawFillingFactor(doFF.codeRawLocalCO12,drawFigure=True,inputID=462793, testPoly=testPoly, individuals=drawIndi )

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






if 0:
    doMain.drawChiSquareTest(trainingRatio=0.8 )
    doMain.drawChiSquareTest(trainingRatio=0.7 )

    sys.exit()




if 0:
    doFF.calCode=doFF.codeRawLocalCO12 #codeRawLocalCO12

    doFF.printCloudNumberTable()
    #doFF.printFillingCat()
    #doFF.printFluxCat()


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







if 0:

    doFF.calCode=doFF.codeRawLocalCO12
    doFF.testEqualCutOff()

    sys.exit()


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

if 0:
    #doMain.getCuoffFlux( doFF.codeRawLocalCO18 )
    doMain.getCutffFF( doFF.codeRawLocalCO12 )

    sys.exit()





if 0:
    #doMain.getCuoffFlux( doFF.codeRawLocalCO18 )
    #doMain.getCutffFF( doFF.codeRawLocalCO13 )

    doFF.calCode=doFF.codeRawLocalCO12
    TB=Table.read("cutoffFlux_pureEdgeInfo_fillingFactor_fluxTB_rawLocalCO12rawLocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit")

    TB.sort("sum")
    TB.reverse()


    newTB=doFF.addMWISPFFerrorCutoff(TB)

    print newTB


    #doFF.calculateFillingFactorCutoff(TB,drawFigure=True)
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
     #doFF.drawCloudSizeChange(drawCode="min")
    #doFF.drawCloudSizeChange(drawCode="max")
    #doFF.drawCloudSizeChange(drawCode="mean")

    doFF.drawCloudNumberChange()

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