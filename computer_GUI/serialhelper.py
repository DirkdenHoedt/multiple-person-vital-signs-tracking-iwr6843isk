import serial
import time

class SerialHelper:
    def __init__(self, serialPort, baudrate) -> None:
        self.serialPort = serial.Serial(serialPort, baudrate)

    def write(self, line):
        self.serialPort.write(bytes(line, 'utf-8'))

    def read(self, bytes = 1):
        return str(self.serialPort.read(bytes), 'utf-8')

    def readline(self):
        return str(self.serialPort.readline(), 'utf-8')

    def readall(self):
        return self.serialPort.read_all()

    def sendConfig(self, config):
        # Emtry the user serial port
        print(str(self.serialPort.read_all(), 'utf-8'))

        for line in config:
            if line == '\n' or line == ' \n' or line.startswith('%'):
                continue

            self.write(line)
            print(line, end='')

            for i in range(3):
                ret = self.readline()
                if 'Done' in ret:
                    print(ret, end='')
                    break
                if 'not recognized as a CLI command' in ret:
                    print(ret, end='')
                    return
                if 'Error' in ret:
                    print(ret, end='')
                    return
            time.sleep(0.2)
