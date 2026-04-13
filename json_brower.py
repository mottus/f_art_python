# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 11:51:22 2026
@author: matti.mottus@gmail.com 
created largely using claude
"""

import tkinter as tk
from tkinter import ttk, filedialog
import json
import os

# ── the function that runs on the loaded data ────────────────────────────────
def process(data):
    """
    Replace the body of this function with whatever you want to do.
    'data' is the parsed JSON (dict, list, etc.).
    Return a string to display in the output panel.
    """
    n = len(data) if isinstance(data, (dict, list)) else 1
    return f"Loaded {type(data).__name__} with {n} top-level items.\n\n{json.dumps(data, indent=2)}"

# ── folder / file helpers ────────────────────────────────────────────────────
def browse_folder():
    folder = filedialog.askdirectory(title="Select folder")
    if folder:
        folder_var.set(folder)
        load_files(folder)

def browse_file():
    path = filedialog.askopenfilename(
        title="Select JSON file",
        filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
        initialdir=folder_var.get(),
    )
    if path:
        folder_var.set(os.path.dirname(path))
        load_files(os.path.dirname(path))
        filename = os.path.basename(path)
        items = listbox.get(0, tk.END)
        if filename in items:
            idx = items.index(filename)
            listbox.selection_clear(0, tk.END)
            listbox.selection_set(idx)
            listbox.see(idx)
        open_file(path)

def load_files(folder):
    listbox.delete(0, tk.END)
    try:
        files = sorted(f for f in os.listdir(folder) if f.endswith(".json"))
    except PermissionError:
        return
    for f in files:
        listbox.insert(tk.END, f)
    status_var.set(f"{len(files)} file(s) in {folder}")

def on_listbox_select(event=None):
    selection = listbox.curselection()
    if not selection:
        return
    filename = listbox.get(selection[0])
    open_file(os.path.join(folder_var.get(), filename))

def open_file(path):
    global current_data
    try:
        with open(path) as f:
            current_data = json.load(f)
        json_text.delete("1.0", tk.END)
        json_text.insert(tk.END, json.dumps(current_data, indent=2))
        output_text.delete("1.0", tk.END)
        status_var.set(f"Loaded: {os.path.basename(path)}")
        run_btn.config(state=tk.NORMAL)
    except Exception as e:
        status_var.set(f"Error loading file: {e}")

# ── run button ───────────────────────────────────────────────────────────────
def on_run():
    if current_data is None:
        return
    try:
        result = process(current_data)
        output_text.delete("1.0", tk.END)
        output_text.insert(tk.END, result)
        status_var.set("Done.")
    except Exception as e:
        output_text.delete("1.0", tk.END)
        output_text.insert(tk.END, f"Error:\n{e}")
        status_var.set("Error — see output panel.")

# ── build UI ─────────────────────────────────────────────────────────────────
current_data = None

root = tk.Tk()
root.title("JSON viewer")
root.geometry("950x640")

# status bar
status_var = tk.StringVar(value="No file loaded.")
ttk.Label(root, textvariable=status_var, anchor=tk.W,
          relief=tk.SUNKEN).pack(side=tk.BOTTOM, fill=tk.X)

# ── top bar: folder path + browse buttons ────────────────────────────────────
top = ttk.Frame(root, padding=(10, 8, 10, 0))
top.pack(fill=tk.X)

folder_var = tk.StringVar(value=os.getcwd())
ttk.Label(top, text="Folder:").pack(side=tk.LEFT)
ttk.Entry(top, textvariable=folder_var, width=55).pack(side=tk.LEFT, padx=(4, 4))
ttk.Button(top, text="Browse folder…",
           command=browse_folder).pack(side=tk.LEFT, padx=(0, 4))
ttk.Button(top, text="Open file…",
           command=browse_file).pack(side=tk.LEFT)
ttk.Button(top, text="Refresh",
           command=lambda: load_files(folder_var.get())).pack(side=tk.LEFT, padx=(4, 0))

# ── main area ────────────────────────────────────────────────────────────────
main = ttk.Frame(root, padding=10)
main.pack(fill=tk.BOTH, expand=True)

# left: file list
left = ttk.Frame(main)
left.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))

ttk.Label(left, text="Files").pack(anchor=tk.W)
listbox = tk.Listbox(left, width=22, activestyle="dotbox")
listbox.pack(fill=tk.Y, expand=True)
listbox.bind("<<ListboxSelect>>", on_listbox_select)

# right: JSON + output
right = ttk.Frame(main)
right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

paned = ttk.PanedWindow(right, orient=tk.VERTICAL)
paned.pack(fill=tk.BOTH, expand=True)

# JSON panel
json_frame = ttk.LabelFrame(paned, text="JSON content", padding=4)
paned.add(json_frame, weight=1)

json_text = tk.Text(json_frame, wrap=tk.NONE, font=("Courier", 11))
json_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
ttk.Scrollbar(json_frame, orient=tk.VERTICAL,
              command=json_text.yview).pack(side=tk.RIGHT, fill=tk.Y)

# run button
btn_bar = ttk.Frame(right)
btn_bar.pack(fill=tk.X, pady=6)
run_btn = ttk.Button(btn_bar, text="▶  Run", command=on_run, state=tk.DISABLED)
run_btn.pack(side=tk.LEFT)
ttk.Label(btn_bar, text="← runs process(data) on the loaded file",
          foreground="gray").pack(side=tk.LEFT, padx=8)

# output panel
out_frame = ttk.LabelFrame(paned, text="Output", padding=4)
paned.add(out_frame, weight=1)

output_text = tk.Text(out_frame, wrap=tk.NONE, font=("Courier", 11))
output_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
ttk.Scrollbar(out_frame, orient=tk.VERTICAL,
              command=output_text.yview).pack(side=tk.RIGHT, fill=tk.Y)

load_files(folder_var.get())
root.mainloop()