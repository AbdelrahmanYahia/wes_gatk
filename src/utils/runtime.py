from datetime import datetime

class RunTime:
    def __init__(self, start_time=datetime.now().strftime("%Y-%m-%d %H:%M:%S")):
        self.start_time = datetime.strptime(start_time, '%Y-%m-%d %H:%M:%S')

    def elapsed(self):
        self.end_time = datetime.now()
        return f'Duration: {self.end_time - self.start_time}'