import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

class PlotHelper:
    def __init__(self) -> None:
        # Init all the graphs we are going to use at the right location
        self.axHeatmap = plt.subplot2grid((4, 8), (0, 0), colspan=2, rowspan=4)
        self.axBreathZone1 = plt.subplot2grid((4, 8), (0, 2), colspan=3)
        self.axHeartZone1 = plt.subplot2grid((4, 8), (1, 2), colspan=3)
        self.axBreathZone2 = plt.subplot2grid((4, 8), (0, 5), colspan=3)
        self.axHeartZone2 = plt.subplot2grid((4, 8), (1, 5), colspan=3)
        self.axBreathZone3 = plt.subplot2grid((4, 8), (2, 2), colspan=3)
        self.axHeartZone3 = plt.subplot2grid((4, 8), (3, 2), colspan=3)
        self.axBreathZone4 = plt.subplot2grid((4, 8), (2, 5), colspan=3)
        self.axHeartZone4 = plt.subplot2grid((4, 8), (3, 5), colspan=3)
        
        # Set proper limits for the graphs up front
        self.axBreathZone1.set_ylim([-1, 1])
        self.axHeartZone1.set_ylim([-1, 1])
        self.axBreathZone2.set_ylim([-1, 1])
        self.axHeartZone2.set_ylim([-1, 1])
        self.axBreathZone3.set_ylim([-1, 1])
        self.axHeartZone3.set_ylim([-1, 1])
        self.axBreathZone4.set_ylim([-1, 1])
        self.axHeartZone4.set_ylim([-1, 1])

        # Init all of the data holders
        self.rangeAzimuth = np.zeros((48, 64))
        
        self.breathZone1 = np.zeros(128)
        self.heartZone1 = np.zeros(128)
        self.heartrateZone1 = 0
        self.breathingrateZone1 = 0

        self.breathZone2 = np.zeros(128)
        self.heartZone2 = np.zeros(128)
        self.heartrateZone2 = 0
        self.breathingrateZone2 = 0

        self.breathZone3 = np.zeros(128)
        self.heartZone3 = np.zeros(128)
        self.heartrateZone3 = 0
        self.breathingrateZone3 = 0

        self.breathZone4 = np.zeros(128)
        self.heartZone4 = np.zeros(128)
        self.heartrateZone4 = 0
        self.breathingrateZone4 = 0

        # Populate the graphs with the right type
        plt.ion()
        self.graphHeatmap = self.axHeatmap.imshow(np.random.rand(64, 48), vmin=0, vmax=1200)
        self.graphBreathZone1, = self.axBreathZone1.plot(self.breathZone1, 'b')
        self.graphHeartZone1, = self.axHeartZone1.plot(self.heartZone1, 'r')
        self.graphBreathZone2, = self.axBreathZone2.plot(self.breathZone2, 'b')
        self.graphHeartZone2, = self.axHeartZone2.plot(self.heartZone2, 'r')
        self.graphBreathZone3, = self.axBreathZone3.plot(self.breathZone3, 'b')
        self.graphHeartZone3, = self.axHeartZone3.plot(self.heartZone3, 'r')
        self.graphBreathZone4, = self.axBreathZone4.plot(self.breathZone4, 'b')
        self.graphHeartZone4, = self.axHeartZone4.plot(self.heartZone4, 'r')

        self.axBreathZone1.set_title('Red zone')
        self.axBreathZone2.set_title('Yellow zone')
        self.axBreathZone3.set_title('Green zone')
        self.axBreathZone4.set_title('Blue zone')
        self.axBreathZone1.xaxis.set_label_position('top') 
        self.axBreathZone2.xaxis.set_label_position('top') 

        self.rects = []
        cols = ['r', 'y', 'g', 'b']
        for i in range(4):
            rect = patches.Rectangle((20, 20), 6, 6, linewidth=1, edgecolor=cols[i], facecolor='none')
            self.axHeatmap.add_patch(rect)
            self.rects.append(rect)

        # Init the framecounter
        self.frameCounter = 0
    
    def update(self, heatmapData, vitalSigns, zones):
        # Update the heatmap
        self.graphHeatmap.set_data(heatmapData)

        # Add the vital sign data to the right dataholder and update graph
        self.breathZone1[self.frameCounter] = vitalSigns[1]
        self.heartZone1[self.frameCounter] = vitalSigns[2]
        self.heartrateZone1 = vitalSigns[3]
        self.breathingrateZone1 = vitalSigns[4]
        self.graphBreathZone1.set_ydata(self.breathZone1)
        self.graphHeartZone1.set_ydata(self.heartZone1)
        self.axBreathZone1.set_xlabel("Heartrate: {}, breathing rate: {}".format(np.floor(self.heartrateZone1), np.floor(self.breathingrateZone1)))

        self.breathZone2[self.frameCounter] = vitalSigns[6]
        self.heartZone2[self.frameCounter] = vitalSigns[7]
        self.heartrateZone2 = vitalSigns[8]
        self.breathingrateZone2 = vitalSigns[9]
        self.graphBreathZone2.set_ydata(self.breathZone2)
        self.graphHeartZone2.set_ydata(self.heartZone2)
        self.axBreathZone2.set_xlabel("Heartrate: {}, breathing rate: {}".format(np.floor(self.heartrateZone2), np.floor(self.breathingrateZone2)))

        self.breathZone3[self.frameCounter] = vitalSigns[11]
        self.heartZone3[self.frameCounter] = vitalSigns[12]
        self.heartrateZone3 = vitalSigns[13]
        self.breathingrateZone3 = vitalSigns[14]
        self.graphBreathZone3.set_ydata(self.breathZone3)
        self.graphHeartZone3.set_ydata(self.heartZone3)
        self.axHeartZone3.set_xlabel("Heartrate: {}, breathing rate: {}".format(np.floor(self.heartrateZone3), np.floor(self.breathingrateZone3)))

        self.breathZone4[self.frameCounter] = vitalSigns[16]
        self.heartZone4[self.frameCounter] = vitalSigns[17]
        self.heartrateZone4 = vitalSigns[18]
        self.breathingrateZone4 = vitalSigns[19]
        self.graphBreathZone4.set_ydata(self.breathZone4)
        self.graphHeartZone4.set_ydata(self.heartZone4)
        self.axHeartZone4.set_xlabel("Heartrate: {}, breathing rate: {}".format(np.floor(self.heartrateZone4), np.floor(self.breathingrateZone4)))

        for i, c in enumerate(zones):
            self.rects[i].set_xy((c[1] - 3, c[0] - 3))

        self.frameCounter += 1
        if self.frameCounter > 127:
            self.frameCounter = 0

        plt.draw()
        plt.pause(0.01)