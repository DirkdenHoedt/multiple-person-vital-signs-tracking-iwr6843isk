import re

def read_config_file(path):
    res = []
    with open(path, 'r') as f:
        for line in f:
            res.append(line)
    return res

def pow2roundup(x):
    y = 1
    while x > y:
        y = y * 2
    return y

def parse_config_file(file):
    res = {}
    for l in file:
        l = re.split(' +', l)
        if l[0] == 'channelCfg':
            tx = int(l[2])
            res['numTxAzimAnt'] = (tx & 0b001) + (tx >> 2 & 0b001)
            res['numTxElevAnt'] = 0
            rx = int(l[1])
            res['numRxAnt'] = (rx & 1) + ((rx >> 1) & 1) + ((rx >> 2) & 1) + ((rx >> 3) & 1)
            res['numTxAnt'] = res['numTxAzimAnt'] + res['numTxElevAnt']
        elif l[0] == 'profileCfg':
            res['startFreq'] = int(l[2])
            res['idleTime'] = int(l[3])
            res['rampEndTime'] = int(l[5])
            res['freqSlopeConst'] = int(l[8])
            res['numAdcSamples'] = int(l[10])
            res['digOutSampleRate'] = int(l[11])
        elif l[0] == 'frameCfg':
            res['chirpStartIdx'] = int(l[1])
            res['chirpEndIdx'] = int(l[2])
            res['numLoops'] = int(l[3])
            res['numFrames'] = int(l[4])
            res['framePeriodicity'] = int(l[5])
        elif l[0] == 'guiMonitor':
            res['decision'] = int(l[1])
            res['rangeAzimuthHeatMap'] = int(l[2])
        elif l[0] == 'zoneDef':
            res['numZones'] = int(l[1])
            zoneDef = []
            for i in range(0, res['numZones'] * 4, 4):
                zoneDef.append([int(l[i+2]), int(l[i+3]), int(l[i+4]), int(l[i+5])])
            res['zoneDef'] = zoneDef
        elif l[0] == 'coeffMatrixRow':
            if res.get('coeffMatrix') is None:
                res['coeffMatrix'] = []
            res['coeffMatrix'].append([float(l[3]), float(l[4]), float(l[5]), float(l[6]), float(l[7]), float(l[8])])
        elif l[0] == 'meanVector':
            res['meanVector'] = [float(l[1]), float(l[2]), float(l[3]), float(l[4]), float(l[5]), float(l[6])]
        elif l[0] == 'stdVector':
            res['stdVector'] = [float(l[1]), float(l[2]), float(l[3]), float(l[4]), float(l[5]), float(l[6])]
        elif l[0] == 'oddemoParms':
            res['windowLen'] = int(l[1])
            res['diagLoadFactor'] = float(l[2])
        elif l[0] == 'rowNoise':
            if res.get('rowNoise') is None:
                res['rowNoise'] = []
            for i in range(3, 4 + int(l[2]) - 1):
                res['rowNoise'].append(float(l[i]))

    res['numChirpsPerFrame'] = (res['chirpEndIdx'] - res['chirpStartIdx'] + 1) * res['numLoops']
    res['numDopplerBins'] = res['numChirpsPerFrame'] / res['numTxAnt']
    res['numRangeBins'] = pow2roundup(res['numAdcSamples'])
    res['rangeResolutionMeters'] = 3e8 * res['digOutSampleRate'] * 1e3 / (2 * res['freqSlopeConst'] * 1e12 * res['numAdcSamples'])
    res['rangeIdxToMeters'] = 3e8 * res['digOutSampleRate'] * 1e3 / (2 * res['freqSlopeConst'] * 1e12 * res['numRangeBins'])
    res['dopplerResolutionMps'] = 3e8 / (2 * res['startFreq'] * 1e9 * (res['idleTime'] + res['rampEndTime']) * 1e-6 * res['numDopplerBins'] * res['numTxAnt'])
    res['numAngleBins'] = 48

    return res