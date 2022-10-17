from datetime import datetime

class LogHelper():
    def __init__(self) -> None:
        self.f = open('logs/vs_log_' + datetime.now().strftime('%Y-%m-%d#%H:%M:%S') + '.txt', 'w')

    def write(self, data):
        res = ''
        for d in data:
            res += str(d) + ','
        res[-1] = '\n'
        self.f.write(res)

    def close(self):
        self.f.close()