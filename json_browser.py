# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 11:51:22 2026
@author: matti.mottus@gmail.com 
created largely using claude
"""
import numpy as np
import tkinter as tk
from tkinter import ttk, filedialog
import json
import os
import io
from pathlib import Path
from f_art.frtclass import frt_model
import matplotlib.pyplot as plt

spectradir = str( Path(r'C:\Users\mmattim\2026\kasikirjad\manuscript_FisherIndex\f_art_python\kaikki_spektrit') )
# the spectral plot to be updated during the script. A new one is created when closed
fig = None
ax = None
colorN = 0


def process(data):
    global fig, ax, colorN
    newplot = False
    filename = os.path.basename(current_path)   # e.g. "a_test.json"
    # folder   = os.path.dirname(current_path)    # e.g. "C:/Users/matti/data"

    F = frt_model( frt_datadir=spectradir, frtconf=data )
    F.reflectance()
    printstring = io.StringIO()
    print(F, file=printstring)
    
    # -- plot R, T
    cmap = plt.get_cmap('tab20') 
    if fig is None or not plt.fignum_exists(fig.number):
        fig, ax = plt.subplots(figsize=(8, 4))
        newplot = True
        colorN = 0
    ax.plot(F.wl, F.R,  label='reflectance '+filename,  color=cmap(colorN), lw=1.5)
    ax.plot(F.wl, F.T, label='transmittance '+filename, color=cmap(colorN), lw=1.5, linestyle='--')
    colorN = colorN+1
    if not newplot:
        all_y = [v for line in ax.get_lines() for v in line.get_ydata()
                 if not np.isnan(v)]
        if all_y:
            ax.set_ylim( 0, max(all_y) )
        fig.canvas.draw()
   
    ax.set_xlabel('wavelength (nm)')
    ax.set_ylabel('reflectance,transmittance')
    # ax.set_xlim(400, 1000)
    ax.set_ylim(bottom=0)
    
    ax.legend(frameon=False)
    # ax.set_title(filename)
    fig.tight_layout()
    plt.show()
    
    return printstring.getvalue()

# ── tree view helpers ─────────────────────────────────────────────────────────
def populate_tree(tree, parent, key, value):
    """Recursively insert a JSON node. Leaf nodes store their data path as tags."""
    if isinstance(value, dict):
        label = f"{key}  {{}}" if key != "" else "{}"
        node = tree.insert(parent, tk.END, text=label, open=False)
        for k, v in value.items():
            populate_tree(tree, node, k, v)
    elif isinstance(value, list):
        label = f"{key}  [{len(value)}]" if key != "" else f"[{len(value)}]"
        node = tree.insert(parent, tk.END, text=label, open=False)
        for i, v in enumerate(value):
            populate_tree(tree, node, str(i), v)
    else:
        label = f"{key}:  {value}"
        tree.insert(parent, tk.END, text=label, tags=("leaf",))
 
 
def load_tree(data):
    tree.delete(*tree.get_children())
    if isinstance(data, dict):
        for k, v in data.items():
            populate_tree(tree, "", k, v)
    elif isinstance(data, list):
        for i, v in enumerate(data):
            populate_tree(tree, "", str(i), v)
    else:
        tree.insert("", tk.END, text=str(data))
 
 
# ── resolve a tree item to its key path in current_data ──────────────────────
def item_path(item):
    """Return list of keys/indices leading to this item in current_data."""
    parts = []
    while item:
        text = tree.item(item, "text")
        parent = tree.parent(item)
        if ":" in text:                         # leaf:  "key:  value"
            key = text.split(":")[0].strip()
        else:                                   # branch: "key  {}" or "key  [N]"
            key = text.split("  ")[0].strip()
        if key:
            parts.append(key)
        item = parent
    parts.reverse()
    return parts
 
 
def get_data_at(path):
    """Walk current_data along path, return (parent_container, final_key)."""
    node = current_data
    for key in path[:-1]:
        node = node[int(key)] if isinstance(node, list) else node[key]
    final = path[-1]
    return node, int(final) if isinstance(node, list) else final
 
 
# ── inline editor ─────────────────────────────────────────────────────────────
edit_entry = None   # the floating Entry widget when editing
 
def start_edit(item):
    global edit_entry
    if edit_entry:
        commit_edit()
 
    text = tree.item(item, "text")
    if ":" not in text:          # branch node — not editable
        return
    key, _, val_str = text.partition(":  ")
 
    bbox = tree.bbox(item)
    if not bbox:
        return
    x, y, w, h = bbox
 
    edit_entry = tk.Entry(tree, font=("Courier", 11))
    edit_entry.place(x=x, y=y, width=max(w, 200), height=h)
    edit_entry.insert(0, val_str)
    edit_entry.select_range(0, tk.END)
    edit_entry.focus_set()
    edit_entry._item = item
 
    edit_entry.bind("<Return>",  lambda e: commit_edit())
    edit_entry.bind("<Escape>",  lambda e: cancel_edit())
    edit_entry.bind("<FocusOut>", lambda e: commit_edit())
 
 
def commit_edit():
    global edit_entry
    if not edit_entry:
        return
    item     = edit_entry._item
    new_str  = edit_entry.get()
    edit_entry.destroy()
    edit_entry = None
 
    # parse new value to the right Python type
    try:
        new_val = json.loads(new_str)   # handles int, float, bool, null, string
    except json.JSONDecodeError:
        new_val = new_str               # keep as plain string
 
    # update current_data
    path = item_path(item)
    if not path:
        return
    try:
        container, key = get_data_at(path)
        container[key] = new_val
    except (KeyError, IndexError, TypeError):
        status_var.set("Could not update value.")
        return
 
    # refresh just this node's label (faster than reloading the whole tree)
    orig_text = tree.item(item, "text")
    label_key = orig_text.split(":")[0]
    tree.item(item, text=f"{label_key}:  {new_val}")
    status_var.set(f"Updated {label_key} → {new_val}")
 
 
def cancel_edit():
    global edit_entry
    if edit_entry:
        edit_entry.destroy()
        edit_entry = None
 
 
def on_tree_double_click(event):
    item = tree.identify_row(event.y)
    if item:
        start_edit(item)
 
 
# ── callbacks ────────────────────────────────────────────────────────────────
def on_close():
    root.quit()
    root.destroy()
 
 
def browse_folder():
    folder = filedialog.askdirectory(
        title="Select folder",
        initialdir=folder_entry.cget("text"),
    )
    if folder:
        folder_entry.config(text=folder)
        load_files(folder)
 
 
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
    open_file(os.path.join(folder_entry.cget("text"), filename))
 
 
def open_file(path):
    global current_data, current_path
    try:
        with open(path) as f:
            current_data = json.load(f)
        load_tree(current_data)
        output_text.delete("1.0", tk.END)
        current_path = path
        status_var.set(f"Loaded: {os.path.basename(path)}")
        run_btn.config(state=tk.NORMAL)
    except Exception as e:
        status_var.set(f"Error: {e}")
 
 
def on_save_as():
    if current_data is None:
        return
    folder    = os.path.dirname(current_path) if current_path else folder_entry.cget("text")
    name      = os.path.basename(current_path) if current_path else "output.json"
    stem, ext = os.path.splitext(name)
    suggested = stem + "_modified" + ext
    path = filedialog.asksaveasfilename(
        title="Save as",
        initialdir=folder,
        initialfile=suggested,
        defaultextension=".json",
        filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
    )
    if not path:
        return
    with open(path, "w") as f:
        json.dump(current_data, f, indent=2)
    status_var.set(f"Saved: {os.path.basename(path)}")
 
 
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
current_path = None
 
root = tk.Tk()
root.title("JSON viewer / editor")
root.geometry("1000x700")
root.protocol("WM_DELETE_WINDOW", on_close)
 
# status bar
status_var = tk.StringVar(value="No file loaded.")
ttk.Label(root, textvariable=status_var, anchor=tk.W,
          relief=tk.SUNKEN).pack(side=tk.BOTTOM, fill=tk.X)
 
# top bar
top = ttk.Frame(root, padding=(10, 8, 10, 0))
top.pack(fill=tk.X)
top.columnconfigure(1, weight=1)   # column 1 (the entry) expands
 
try:
    _default_folder = os.path.dirname(os.path.abspath(__file__))
except NameError:
    _default_folder = os.getcwd()
 
ttk.Label(top, text="Folder:").grid(row=0, column=0, sticky=tk.W, padx=(0, 4))
folder_entry = ttk.Label(top, text="", relief="sunken", anchor=tk.W)
folder_entry.grid(row=0, column=1, sticky=tk.EW, padx=(0, 4))
ttk.Button(top, text="Browse…", command=browse_folder).grid(row=0, column=2, padx=(0, 4))
ttk.Button(top, text="Refresh",
           command=lambda: load_files(folder_entry.cget("text"))).grid(row=0, column=3)
 
# main area
main = ttk.Frame(root, padding=10)
main.pack(fill=tk.BOTH, expand=True)
 
# left: file list
left = ttk.Frame(main)
left.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
 
ttk.Label(left, text="Files").pack(anchor=tk.W)
listbox = tk.Listbox(left, width=22, activestyle="dotbox")
listbox.pack(fill=tk.Y, expand=True)
listbox.bind("<<ListboxSelect>>", on_listbox_select)
 
# right column
right = ttk.Frame(main)
right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
 
# button bar — packed first so it is never squeezed out
btn_bar = ttk.Frame(right)
btn_bar.pack(side=tk.BOTTOM, fill=tk.X, pady=6)
run_btn = ttk.Button(btn_bar, text="▶  Run FRT", command=on_run, state=tk.DISABLED)
run_btn.pack(side=tk.LEFT, padx=(0, 6))
ttk.Button(btn_bar, text="Save as…", command=on_save_as).pack(side=tk.LEFT, padx=(0, 6))
ttk.Button(btn_bar, text="Close", command=on_close).pack(side=tk.LEFT)
 
# paned window
paned = ttk.PanedWindow(right, orient=tk.VERTICAL)
paned.pack(fill=tk.BOTH, expand=True)
 
# JSON tree panel
tree_frame = ttk.LabelFrame(paned, text="JSON content  (double-click a value to edit)", padding=4)
paned.add(tree_frame, weight=2)
 
tree = ttk.Treeview(tree_frame, show="tree", selectmode="browse")
tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
ttk.Scrollbar(tree_frame, orient=tk.VERTICAL,
              command=tree.yview).pack(side=tk.RIGHT, fill=tk.Y)
 
tree.bind("<Double-1>", on_tree_double_click)
 
# right-click context menu for TreeClasses items
tree_menu = tk.Menu(root, tearoff=0)
tree_menu.add_command(label="Delete this tree class",
                      command=lambda: delete_tree_class())
 
 
def get_tree_class_index(item):
    parent = tree.parent(item)
    if not parent:
        return None
    if not tree.item(parent, "text").startswith("TreeClasses"):
        return None
    return list(tree.get_children(parent)).index(item)
 
 
def delete_tree_class():
    sel = tree.selection()
    if not sel:
        return
    idx = get_tree_class_index(sel[0])
    if idx is None:
        return
    if messagebox.askyesno("Delete", f"Delete tree class {idx}?"):
        current_data["TreeClasses"].pop(idx)
        load_tree(current_data)
        status_var.set(f"Deleted tree class {idx}.")
 
 
def on_tree_right_click(event):
    item = tree.identify_row(event.y)
    if not item:
        return
    tree.selection_set(item)
    if get_tree_class_index(item) is not None:
        tree_menu.tk_popup(event.x_root, event.y_root)
 
 
tree.bind("<Button-3>", on_tree_right_click)
 
# output panel
out_frame = ttk.LabelFrame(paned, text="Output", padding=4)
paned.add(out_frame, weight=1)
 
output_text = tk.Text(out_frame, wrap=tk.NONE, font=("Courier", 11))
output_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
ttk.Scrollbar(out_frame, orient=tk.VERTICAL,
              command=output_text.yview).pack(side=tk.RIGHT, fill=tk.Y)
 
load_files(folder_entry.cget("text"))
root.mainloop()