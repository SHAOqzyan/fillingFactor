from fillingfactor import checkFillingFactor
import os

doFF=checkFillingFactor()





if 1:#pipeline line

    doFF.calCode="Q1Sgr" #30-70 km/s
    doFF.rawFITS="./data/G2650Sgr70.fits"

    #step 1 smooth by factor
    doFF.smoothFITSbySMFactor(doFF.rawFITS)

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