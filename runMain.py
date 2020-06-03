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


    def getFillingFactorAndEdge(self,calCode):

        #get fluxTB
        doFF.calCode=calCode
        smTBFile = doFF.getSmoothAndNoiseFITSSingle(getCleanTBFile=True)

        fluxTBFile="fluxTB_"+os.path.basename(smTBFile)

        if not os.path.isfile( fluxTBFile):
            print "flux TB does not exist, please check your data or calcode"

            return


        #getFillingfactor

        ffTBName=doFF.calculateFillingFactor(fluxTBFile, drawFigure=False)

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

        return "edgeInfo_fillingFactor_fluxTB_"+ os.path.basename( smTBFile )


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





    def zzz(self):
        """

        :return:
        """
        pass





doMain=fillingMain()

if 1: #check particular clouds, smaller filling factors

    checkTB="edgeInfo_fillingFactor_fluxTB_rawLocalCO12rawLocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"

    testCode=doFF.codeRawLocalCO12

    getTB= doFF.selectBySizeAndFactorRange(checkTB,sizeRange=[0, 4 ],factorRange=[0.6,1])
    getTB.sort(doFF.ffMWISPCol)
    print getTB


    #testIDList=[2607,10959,6656,117419,29631]
    testIDList=[274907, 215799]

    for eachID in testIDList:

        doFF.calFFByID(testCode,eachID,drawFigure=True,useSigmaCut=True)
        doFF.drawInMap(testCode,eachID)


    sys.exit()




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




if 0:


    for eachCode in doFF.allRawCodeList:


        doMain.getFillingFactorAndEdge(eachCode)

    sys.exit()


if 0:  # find clouds #


    for eachCode in doFF.allRawCodeList:
        doFF.getFluxListForEachCloud(calCode= eachCode )

    sys.exit()




if 0: #clean the first fits, to do the flux calculation

    for eachCode in doFF.allRawCodeList:

        doMain.cleanFITS( eachCode, onlyFirst=True )

if 0: # clean all fits, #Do this tonight

    for eachCode in doFF.allRawCodeList:
        doMain.cleanFITS( eachCode, onlyFirst=False )
        doMain.removeUselessFITS(eachCode)

    sys.exit()

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
    doFF.testEqualCutOff()


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