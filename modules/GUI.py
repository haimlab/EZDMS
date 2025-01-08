import tkinter as tk

def on_button_click():
    label.config(text="Hello, World!")

root = tk.Tk()
root.title("Simple Tkinter App")

label = tk.Label(root, text="Press the button!")
label.pack(pady=20)

button = tk.Button(root, text="Click Me", command=on_button_click)
button.pack(pady=20)

root.mainloop()
