import struct
import numpy as np

def getHeader(data, index):
    header = {}
    index += 8
    header['totalPacketLen'] = int.from_bytes(data[index:index+4], 'little')
    index += 4
    header['platform'] = int.from_bytes(data[index:index+4], 'little')
    index += 4
    header['frameNumber'] = int.from_bytes(data[index:index+4], 'little')
    index += 4
    header['timeCpuCycles'] = int.from_bytes(data[index:index+4], 'little')
    index += 4
    header['numDetectedObj'] = int.from_bytes(data[index:index+4], 'little')
    index += 4
    header['numTLVs'] = int.from_bytes(data[index:index+4], 'little')
    index += 4
    return header, index

def getTlv(data, index):
    tlv = {}
    tlv['type'] = int.from_bytes(data[index:index+4], 'little')
    index += 4
    tlv['length'] = int.from_bytes(data[index:index+4], 'little')
    index += 4
    return tlv, index

def getOccupDemoRangeAzimuthHeatMap(data, index, numRangeBins, numAngleBins):
    dataLength = numRangeBins * numAngleBins
    rangeAzimuth = struct.unpack('<' + str(dataLength) + 'I', data[index:index+dataLength*4])
    return np.array(rangeAzimuth).reshape((numRangeBins, numAngleBins)), index + dataLength * 4

def getOccupDemoShortHeatMap(data, index, numRangeBins, numAngleBins):
    dataLength = numRangeBins * numAngleBins
    # print(index, index + dataLength*2)
    rangeAzimuth = struct.unpack('<' + str(dataLength) + 'H', data[index:index+dataLength*2])
    return np.array(rangeAzimuth).reshape((numRangeBins, numAngleBins)), index + dataLength * 2

def getOccupDemoByteHeatMap(data, index, numRangeBins, numAngleBins):
    dataLength = numRangeBins * numAngleBins
    rangeAzimuth = struct.unpack('<' + str(dataLength) + 'B', data[index:index+dataLength])
    return np.array(rangeAzimuth).reshape((numRangeBins, numAngleBins)), index + dataLength

def getOccupDemoDecision(data, index, numZones):
    decisionValue = data[index:index+numZones]
    index += numZones
    return decisionValue, index

def getVitalSignsDemoHeartBreathingRate(data, index):
    # 10 fields and 4 bytes per field
    dataLength = 20 * 4
    values = struct.unpack('20f', data[index:index+dataLength])
    return np.array(values), index + dataLength

def getUpdatedZones(data, index):
    dataLength = 9
    values = struct.unpack('<9B', data[index:index+dataLength])
    return values[8], np.array(values[0:8]).reshape(4, 2), index + dataLength

def dumpRowNoiseValues(data, index, numRangeBins):
    return index + numRangeBins * 4

def getStatsInfo(data, index):
    statsInfo = {}
    statsInfo['interFrameProcessingTime'] = int.from_bytes(data[index:index+4], 'little')
    index += 4
    statsInfo['transmitOutputTime'] = int.from_bytes(data[index:index+4], 'little')
    index += 4
    statsInfo['interFrameProcessingMargin'] = int.from_bytes(data[index:index+4], 'little')
    index += 4
    statsInfo['interChirpProcessingMargin'] = int.from_bytes(data[index:index+4], 'little')
    index += 4
    statsInfo['activeFrameCPULoad'] = int.from_bytes(data[index:index+4], 'little')
    index += 4
    statsInfo['interFrameCPULoad'] = int.from_bytes(data[index:index+4], 'little')
    index += 4
    return statsInfo, index