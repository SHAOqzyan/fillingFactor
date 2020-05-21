from fillingfactor import checkFillingFactor
import os
import sys


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


    def drawFillingFactor(self,calCode,targetNoise):
        """

        :return:
        """
        doFF.calCode = calCode

        doFF.drawFillingFactor(absK=targetNoise,useArea=False )




doMain=fillingMain()

if 1:
    doFF.getFFForEachCloud(doFF.codeLocalCO12, drawFigure=True, useSigmaCut=True)
    #doFF.getFillingFactorByCloudID(doFF.codeLocalCO12,8570)
    sys.exit()



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
    doFF.printFluxCat()
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


if 0: #out CO 12

    #doFF.calCode= doFF.codeOutCO12
    #doFF.smoothFITSbySMFactor(doFF.outCO12FITS)
    #sys.exit()

    doFF.calCode= doFF.codeOutCO13
    doFF.smoothFITSbySMFactor(doFF.outCO13FITS)


    doFF.calCode= doFF.codeOutCO18
    doFF.smoothFITSbySMFactor(doFF.outCO18FITS)


    ##############Sgr
    doFF.calCode= doFF.codeSgrCO12
    doFF.smoothFITSbySMFactor(doFF.sgrCO12FITS)


    doFF.calCode= doFF.codeSgrCO13
    doFF.smoothFITSbySMFactor(doFF.sgrCO13FITS)

    doFF.calCode= doFF.codeSgrCO18
    doFF.smoothFITSbySMFactor(doFF.sgrCO18FITS)

    sys.exit()


if 0:

    ##### Scu
    doFF.calCode= doFF.codeScuCO12
    doFF.smoothFITSbySMFactor(doFF.scuCO12FITS)

    doFF.calCode= doFF.codeScuCO13
    doFF.smoothFITSbySMFactor(doFF.scuCO13FITS)

    doFF.calCode = doFF.codeScuCO18
    doFF.smoothFITSbySMFactor(doFF.scuCO18FITS)

    ##### Local
    doFF.calCode = doFF.codeLocalCO18
    doFF.smoothFITSbySMFactor(doFF.localCO18FITS)

    doFF.calCode= doFF.codeLocalCO12
    doFF.smoothFITSbySMFactor(doFF.localCO12FITS)

    doFF.calCode= doFF.codeLocalCO13
    doFF.smoothFITSbySMFactor(doFF.localCO13FITS)

    sys.exit()



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
    doFF.calCode="Q1Local"

    testFile=   "/home/qzyan/WORK/myDownloads/fillingFactor/tmpFiles/Q1LocalG2650Local30_SmFactor_2.0_noise_0.5absK.fits"
    testTBFile= "/home/qzyan/WORK/myDownloads/fillingFactor/tmpFiles/Q1LocalG2650Local30_SmFactor_2.0_noise_0.5absKdbscanS2P4Con1_Clean.fit"

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