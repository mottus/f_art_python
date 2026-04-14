# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 11:51:22 2026
@author: matti.mottus@gmail.com 
created largely using claude
"""
import numpy as np
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import json
import os
import io
from pathlib import Path
from f_art.frtclass import frt_model
import matplotlib.pyplot as plt

current_dir = Path(__file__).resolve().parent

# MODIFY THIS AS NEEDED:
spectradir = str( current_dir / "data" )

# variables for conistent plotting
# the spectral plot to be updated during the script. A new one is created when closed
fig = None
ax = None
colorN = 0


# ── runs when the Run FRT button is pressed ───────────────────────────────────

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

# Pure Claude code from here onwards

# ── runs automatically every time a json file is loaded ──────────────────────
def on_load(data):
    """
    Prints a summary table of tree classes.
    """
    import math
 
    tree_classes = data.get("TreeClasses", [])
    if not tree_classes:
        return "No TreeClasses found."
 
    W       = [28, 12, 10, 10, 16]
    names   = ["Description",  "Density",  "DBH",    "Height",    "Basal area"]
    units   = ["",             "(m⁻²)",    "(cm)",   "(m)",       "(m²/ha)"]
    sep     = "  ".join("-" * w for w in W)
    fmt_row = lambda vals: "  ".join(f"{v:<{w}}" for v, w in zip(vals, W))
 
    lines = [fmt_row(names), fmt_row(units), sep]
    for i, tc in enumerate(tree_classes):
        desc    = str(tc.get("Description", i + 1))
        density = tc.get("StandDensity", None)
        dbh     = tc.get("DBH", None)
        height  = tc.get("TreeHeight", None)
 
        if density is not None and dbh is not None:
            ba     = density * math.pi * (dbh / 200.0) ** 2 * 10000
            ba_str = f"{ba:.1f}"
        else:
            ba_str = "N/A"
 
        density_str = f"{density:.4f}" if density is not None else "N/A"
        dbh_str     = f"{dbh:.1f}"     if dbh     is not None else "N/A"
        height_str  = f"{height:.1f}"  if height  is not None else "N/A"
 
        lines.append(fmt_row([desc, density_str, dbh_str, height_str, ba_str]))
 
    return "\n".join(lines)
 
 
# ── edit-state helpers ───────────────────────────────────────────────────────
_JSON_FRAME_BASE = "JSON content  (double-click a value to edit)"
 
def mark_edited():
    """
    Called whenever current_data is modified interactively.
    - Adds '  [edited]' to the JSON content frame title so the user knows
      there are unsaved changes.
    - Restores the listbox highlight that gets cleared when focus moves to
      the inline entry widget.
    """
    tree_frame.config(text=_JSON_FRAME_BASE + "  [edited]")
    # restore listbox highlight: find which entry matches the current file
    # and re-select it without triggering on_listbox_select (which would reload)
    if current_path:
        fname = os.path.basename(current_path)
        items = listbox.get(0, tk.END)
        if fname in items:
            idx = list(items).index(fname)
            listbox.selection_clear(0, tk.END)
            listbox.selection_set(idx)
            listbox.see(idx)
 
 
def clear_edited():
    """
    Called when a fresh file is loaded.
    Removes the '[edited]' marker from the frame title.
    """
    tree_frame.config(text=_JSON_FRAME_BASE)
 
 
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
 
 
def is_list_node(item):
    """Return the list's key name if this tree item is a list branch, else None."""
    text = tree.item(item, "text")
    # list branch labels look like  "key  [N]"
    if "  [" in text and text.endswith("]") and ":" not in text:
        key = text.split("  [")[0].strip()
        return key if key else None
    return None
 
 
def replace_list_with_file(item):
    """
    Called when the user selects "Replace list with file name…" from the
    right-click menu.
 
    Behaviour:
      - The original list element (e.g. "wl") is DELETED from the container.
      - A new element named "wlFile" (= original key + "File") is ADDED,
        holding the filename string the user types.
      - The tree is refreshed to reflect the change.
 
    An inline Entry widget is placed over the tree node so the user can
    type the filename and confirm with Enter (or cancel with Escape).
    """
    # determine the key name of the list node ("wl", "SQratio", etc.)
    key = is_list_node(item)
    if key is None:
        return   # not a list node — nothing to do
 
    # resolve the path through current_data to find the parent container
    path = item_path(item)
    if not path:
        return
    try:
        container, data_key = get_data_at(path)
    except (KeyError, IndexError, TypeError):
        return
 
    # the new element name: "wl" → "wlFile"
    file_key  = key + "File"
    suggested = file_key + ".txt"   # pre-fill with a sensible default
 
    # place an inline Entry widget over the tree node
    bbox = tree.bbox(item)
    if not bbox:
        return
    x, y, w, h = bbox
 
    global edit_entry
    if edit_entry:
        commit_edit()   # commit any other open edit first
 
    edit_entry = tk.Entry(tree, font=("Courier", 11))
    edit_entry.place(x=x, y=y, width=max(w, 300), height=h)
    edit_entry.insert(0, suggested)
    edit_entry.select_range(0, tk.END)
    edit_entry.focus_set()
 
    # store everything the commit callback will need on the entry widget itself
    edit_entry._item      = item
    edit_entry._container = container   # dict that owns the list
    edit_entry._list_key  = data_key    # original key ("wl")
    edit_entry._file_key  = file_key    # new key ("wlFile")
 
    def _commit_file(event=None):
        """
        Confirm the replacement:
          1. Delete the original list from the container.
          2. Add a new element <orig_key>File with the filename string.
          3. Reload the tree.
        """
        global edit_entry
        if not edit_entry:
            return   # already committed or cancelled
 
        fname    = edit_entry.get().strip()
        cont     = edit_entry._container
        orig_key = edit_entry._list_key
        fkey     = edit_entry._file_key
 
        edit_entry.destroy()
        edit_entry = None
 
        if not fname:
            return   # empty input → cancel silently
 
        # delete the list element
        del cont[orig_key]
 
        # add the new filename element
        cont[fkey] = fname
 
        # reload the whole tree so the change is visible
        load_tree(current_data)
        status_var.set(f"'{orig_key}' removed, '{fkey}' = '{fname}'")
        mark_edited()
 
    edit_entry.bind("<Return>",   _commit_file)
    edit_entry.bind("<Escape>",   lambda e: cancel_edit())
    edit_entry.bind("<FocusOut>", _commit_file)
 
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
    mark_edited()   # flag unsaved changes and restore listbox highlight
 
 
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
    if not folder or not os.path.isdir(folder):
        return
    listbox.delete(0, tk.END)
    try:
        files = sorted(f for f in os.listdir(folder) if f.endswith(".json"))
    except (PermissionError, FileNotFoundError):
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
        clear_edited()          # remove any [edited] marker from previous file
        output_text.delete("1.0", tk.END)
        current_path = path
        folder_entry.config(text=os.path.dirname(path))
        status_var.set(f"Loaded: {os.path.basename(path)}")
        run_btn.config(state=tk.NORMAL)
        try:
            result = on_load(current_data)
            if result:
                output_text.delete("1.0", tk.END)
                output_text.insert(tk.END, result)
        except Exception as ex:
            output_text.delete("1.0", tk.END)
            output_text.insert(tk.END, f"on_load error:\n{ex}")
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
paned.add(tree_frame, weight=3)
 
tree = ttk.Treeview(tree_frame, show="tree", selectmode="browse")
tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
ttk.Scrollbar(tree_frame, orient=tk.VERTICAL,
              command=tree.yview).pack(side=tk.RIGHT, fill=tk.Y)
 
# Double-clicking a leaf node opens an inline editor.
# <Double-1> is tkinter's name for a left-button double-click.
tree.bind("<Double-1>", on_tree_double_click)
 
 
def get_tree_class_child(item):
    """
    Walk up the tree from 'item' until a direct child of the TreeClasses
    node is found.  Returns that child item, or None if the item is not
    inside TreeClasses at any depth.
 
    This means right-clicking either the "0 {}" header node OR any leaf
    inside it (e.g. "TreeClumping: 1.12") will correctly identify which
    tree class is being targeted.
    """
    current = item
    while current:
        parent = tree.parent(current)
        if not parent:
            return None   # reached the root without finding TreeClasses
        if tree.item(parent, "text").startswith("TreeClasses"):
            return current   # current is "0 {}", "1 {}", etc.
        current = parent
    return None
 
 
def _tc_index(item):
    """
    Return (child_node, index) for the tree class containing 'item',
    where child_node is the direct child of TreeClasses ("0 {}", etc.)
    and index is its 0-based position in current_data["TreeClasses"].
    Returns (None, None) if the item is not inside TreeClasses.
    """
    child = get_tree_class_child(item)
    if child is None:
        return None, None
    siblings = list(tree.get_children(tree.parent(child)))
    return child, siblings.index(child)
 
 
def _tc_has_description(idx):
    """Return True if tree class at index idx already has a Description key."""
    return "Description" in current_data["TreeClasses"][idx]
 
 
def delete_tree_class():
    """
    Delete the tree class that was right-clicked.
 
    tree_menu._target_item holds the item that was under the cursor at
    right-click time (stored there by on_tree_right_click so it is still
    correct when the user eventually clicks the menu entry).
    """
    item = tree_menu._target_item
    if item is None:
        return
    child, idx = _tc_index(item)
    if idx is None:
        return
    current_data["TreeClasses"].pop(idx)   # remove from data
    load_tree(current_data)                # refresh tree view
    status_var.set(f"Deleted tree class {idx}.")
    mark_edited()
 
 
def duplicate_tree_class():
    """
    Append a deep copy of the right-clicked tree class to the end of
    TreeClasses.  The copy is independent — editing it does not affect
    the original.
    """
    import copy
    item = tree_menu._target_item
    if item is None:
        return
    child, idx = _tc_index(item)
    if idx is None:
        return
    original = current_data["TreeClasses"][idx]
    duplicate = copy.deepcopy(original)   # deep copy so nested lists are independent
 
    # if Description exists, append " (copy)" to make it distinguishable
    if "Description" in duplicate:
        duplicate["Description"] = str(duplicate["Description"]) + " (copy)"
 
    current_data["TreeClasses"].append(duplicate)
    load_tree(current_data)
    status_var.set(f"Duplicated tree class {idx} → appended as class {len(current_data['TreeClasses'])-1}.")
    mark_edited()
 
 
def add_description():
    """
    Add a Description key to the right-clicked tree class.
 
    An inline Entry is placed over the tree class node so the user can
    type the description and confirm with Enter.  Only shown when the
    tree class does not already have a Description.
    """
    item = tree_menu._target_item
    if item is None:
        return
    child, idx = _tc_index(item)
    if idx is None:
        return
 
    # place an inline entry over the tree class header node
    bbox = tree.bbox(child)   # use child (the "0 {}" node) for correct position
    if not bbox:
        return
    x, y, w, h = bbox
 
    global edit_entry
    if edit_entry:
        commit_edit()
 
    edit_entry = tk.Entry(tree, font=("Courier", 11))
    edit_entry.place(x=x, y=y, width=max(w, 250), height=h)
    edit_entry.insert(0, f"Class {idx + 1}")   # sensible default
    edit_entry.select_range(0, tk.END)
    edit_entry.focus_set()
    edit_entry._item = child
    edit_entry._tc_idx = idx   # remember which tree class to update
 
    def _commit_desc(event=None):
        global edit_entry
        if not edit_entry:
            return
        desc    = edit_entry.get().strip()
        tc_idx  = edit_entry._tc_idx
        edit_entry.destroy()
        edit_entry = None
        if not desc:
            return   # empty → cancel
        # insert Description as the first key so it appears at the top
        tc = current_data["TreeClasses"][tc_idx]
        updated = {"Description": desc}
        updated.update(tc)               # merge remaining keys after Description
        current_data["TreeClasses"][tc_idx] = updated
        load_tree(current_data)
        status_var.set(f"Added Description '{desc}' to tree class {tc_idx}.")
        mark_edited()
 
    edit_entry.bind("<Return>",   _commit_desc)
    edit_entry.bind("<Escape>",   lambda e: cancel_edit())
    edit_entry.bind("<FocusOut>", _commit_desc)
 
 
# ── right-click context menu ──────────────────────────────────────────────────
# The menu is built once and reused.  The item that was under the cursor is
# stored on the menu object itself (tree_menu._target_item) at right-click
# time, so it is available when the user eventually clicks a menu entry.
# This avoids closure/timing problems that occur if the selection changes
# between the right-click and the menu command firing.
#
# Menu entries and when they are enabled:
#   "Delete this tree class"       → always enabled inside TreeClasses
#   "Duplicate this tree class"    → always enabled inside TreeClasses
#   "Add Description"              → enabled only if Description is missing
#   "Replace list with file name…" → enabled only for list [N] nodes
tree_menu = tk.Menu(root, tearoff=0)
tree_menu._target_item = None
tree_menu.add_command(label="Delete this tree class",    command=delete_tree_class)
tree_menu.add_command(label="Duplicate this tree class", command=duplicate_tree_class)
tree_menu.add_command(label="Add Description",           command=add_description)
tree_menu.add_separator()
tree_menu.add_command(label="Replace list with file name…",
                      command=lambda: replace_list_with_file(tree_menu._target_item))
 
 
def on_tree_right_click(event):
    """
    Called when the user right-clicks anywhere on the tree.
 
    Steps:
      1. Identify which row the cursor is over.
      2. Highlight that row.
      3. Check whether the row is inside TreeClasses (→ tree class options)
         or is a list node (→ Replace list option).
      4. Enable/disable each menu entry accordingly.
      5. Store the item on the menu object so commands can access it later.
      6. Pop up the menu at the cursor position.
 
    The menu is only shown if at least one entry applies.
    Right-clicking on ordinary leaf values or non-list branch nodes shows nothing.
    """
    item = tree.identify_row(event.y)   # row under the cursor (empty if none)
    if not item:
        return   # clicked on empty space — do nothing
 
    tree.selection_set(item)   # highlight the clicked row
 
    # determine what the clicked item is
    child, idx  = _tc_index(item)
    is_tc       = idx is not None                          # inside TreeClasses?
    is_list     = is_list_node(item) is not None           # a list node [N]?
    has_desc    = is_tc and _tc_has_description(idx)       # already has Description?
 
    if is_tc or is_list:
        # save the item now — commands run later when the user clicks the entry
        tree_menu._target_item = item
 
        # enable/disable each entry based on context
        tree_menu.entryconfig("Delete this tree class",
                              state=tk.NORMAL if is_tc   else tk.DISABLED)
        tree_menu.entryconfig("Duplicate this tree class",
                              state=tk.NORMAL if is_tc   else tk.DISABLED)
        tree_menu.entryconfig("Add Description",
                              state=tk.NORMAL if (is_tc and not has_desc) else tk.DISABLED)
        tree_menu.entryconfig("Replace list with file name…",
                              state=tk.NORMAL if is_list else tk.DISABLED)
 
        # show the menu at the screen coordinates of the mouse cursor
        tree_menu.tk_popup(event.x_root, event.y_root)
 
 
# <Button-3> is tkinter's name for a right mouse button click.
# Binding it to on_tree_right_click activates the context menu whenever
# the user right-clicks anywhere on the tree widget.
tree.bind("<Button-3>", on_tree_right_click)
 
# output panel
out_frame = ttk.LabelFrame(paned, text="Output", padding=4)
paned.add(out_frame, weight=2)
 
output_text = tk.Text(out_frame, wrap=tk.NONE, font=("Courier", 9))
output_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
ttk.Scrollbar(out_frame, orient=tk.VERTICAL,
              command=output_text.yview).pack(side=tk.RIGHT, fill=tk.Y)
 
# auto-load json files from the script's own folder on startup
try:
    _startup_folder = os.path.dirname(os.path.abspath(__file__))
except NameError:
    _startup_folder = os.getcwd()
 
if os.path.isdir(_startup_folder):
    folder_entry.config(text=_startup_folder)
    load_files(_startup_folder)
 
def _set_sash():
    root.update_idletasks()
    total = paned.winfo_height()
    if total > 10:
        paned.sashpos(0, int(total * 0.60))
    else:
        root.after(100, _set_sash)   # retry if not drawn yet
 
root.after(200, _set_sash)
root.mainloop()