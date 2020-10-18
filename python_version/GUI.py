import tkinter as tk
from PIL import Image, ImageTk

import multiprocessing

class App:
    def __init__(self, master=tk.Tk()):
        self.master = master
        self.fig_size = [1000, 760]
        self.frame = tk.Frame(master)
        self.canvas = tk.Canvas(self.frame, width=1280, height=800)
        self.canvas.pack()

        self.load_image('test_data/sample_8.jpg')
        self.image_label = tk.Label(self.canvas, image=self.fig_image)
        self.image_label.pack()

        self.button_left = tk.Button(self.frame, text="BUTTON",
                                     command=self.update)
        self.button_left.pack(side="left")

        self.frame.bind("q", self.close)
        self.frame.bind("<Escape>", self.close)
        self.frame.pack()
        self.frame.focus_set()

        self.is_active = True

    def load_image(self, filename):
        self.fig_image = ImageTk.PhotoImage(Image.open(filename).resize(self.fig_size, Image.BILINEAR))

    def update(self, *args):
        self.load_image('test_data/sample_8.jpg')
        self.image_label.config(image=self.fig_image)

    def close(self, *args):
        print('GUI closed...')
        self.master.quit()
        self.is_active = False

    def is_closed(self):
        return not self.is_active

    def mainloop(self):
        self.master.mainloop()
        print('mainloop closed...')

if __name__ == '__main__':
    import time
    app = App()
    app.mainloop()

 