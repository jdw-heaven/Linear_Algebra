import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel
from PyQt5.QtCore import Qt
import subprocess

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.label = QLabel("请划词", self)
        self.label.setAlignment(Qt.AlignCenter)
        self.setCentralWidget(self.label)

        self.setFixedSize(300, 200)
        self.setWindowTitle("划词翻译应用")

    def mousePressEvent(self, event):
        # 执行划词翻译服务的curl命令
        command = 'curl "http://127.0.0.1:60827/selection_translate"'
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        output, _ = process.communicate()

        # 更新界面上的翻译结果
        if output:
            translation = output.decode().strip()
            self.label.setText(translation)
        else:
            self.label.setText("请求失败或无翻译结果")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())