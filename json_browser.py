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

# ── copy/paste clipboard ──────────────────────────────────────────────────────
# clipboard holds: {'key': str, 'value': any, 'is_tree_class': bool}
# key is the top-level JSON key (e.g. 'thetv', 'TreeClasses')
# value is a deepcopy of the data at that key (or the single tree-class dict)
_clipboard = None   # None means empty
_has_unsaved_edits = False   # True when current_data differs from the file on disk


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
    global _has_unsaved_edits
    _has_unsaved_edits = True
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
    global _has_unsaved_edits
    _has_unsaved_edits = False
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
    """Return the list's key name if this tree item is a replaceable list branch, else None.
    Only matches lists whose children are scalar leaves (i.e. data arrays),
    not lists of dicts/sub-lists such as TreeClasses."""
    text = tree.item(item, "text")
    # list branch labels look like  "key  [N]"
    if not ("  [" in text and text.endswith("]") and ":" not in text):
        return None
    key = text.split("  [")[0].strip()
    if not key:
        return None
    # exclude lists whose children are objects or sub-lists (e.g. TreeClasses)
    children = tree.get_children(item)
    if children:
        first_text = tree.item(children[0], "text")
        # object children look like "0  {}" and list children like "0  [N]"
        if first_text.endswith("}") or (first_text.endswith("]") and ":" not in first_text):
            return None
    return key
 
 
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
 
 
def _is_edited():
    """Return True if the current file has unsaved changes."""
    return _has_unsaved_edits


def _confirm_switch(new_path):
    """
    If the current file has unsaved changes, ask the user what to do.
    Returns True if it is safe to proceed with loading new_path, False to abort.

    Three choices:
      Save & switch  — write current_data to current_path in-place, then proceed.
      Discard        — proceed without saving.
      Cancel         — stay on the current file.
    """
    if not _is_edited() or current_data is None or current_path is None:
        return True   # nothing to worry about

    fname = os.path.basename(current_path)
    dlg   = tk.Toplevel(root)
    dlg.title("Unsaved changes")
    dlg.resizable(False, False)
    dlg.grab_set()          # modal
    dlg.focus_set()

    ttk.Label(dlg, text=f'"{fname}" has unsaved changes.',
              padding=(16, 12, 16, 4)).pack()
    ttk.Label(dlg, text="What would you like to do?",
              padding=(16, 0, 16, 12)).pack()

    btn_frame = ttk.Frame(dlg, padding=(12, 0, 12, 14))
    btn_frame.pack()

    choice = tk.StringVar(value="cancel")

    def _pick(val):
        choice.set(val)
        dlg.destroy()

    ttk.Button(btn_frame, text="Save & switch",
               command=lambda: _pick("save")).pack(side=tk.LEFT, padx=(0, 6))
    ttk.Button(btn_frame, text="Save copy & switch",
               command=lambda: _pick("savecopy")).pack(side=tk.LEFT, padx=(0, 6))
    ttk.Button(btn_frame, text="Discard",
               command=lambda: _pick("discard")).pack(side=tk.LEFT, padx=(0, 6))
    ttk.Button(btn_frame, text="Cancel",
               command=lambda: _pick("cancel")).pack(side=tk.LEFT)

    # centre dialog over root
    root.update_idletasks()
    rx, ry = root.winfo_rootx(), root.winfo_rooty()
    rw, rh = root.winfo_width(), root.winfo_height()
    dlg.update_idletasks()
    dw, dh = dlg.winfo_width(), dlg.winfo_height()
    dlg.geometry(f"+{rx + (rw - dw)//2}+{ry + (rh - dh)//2}")

    dlg.wait_window()

    if choice.get() == "save":
        try:
            with open(current_path, "w") as f:
                json.dump(current_data, f, indent=2)
            status_var.set(f"Saved: {os.path.basename(current_path)}")
        except OSError as e:
            status_var.set(f"Save error: {e}")
            return False
        return True
    elif choice.get() == "savecopy":
        new_fname = _write_copy(current_path, current_data)
        if new_fname:
            load_files(os.path.dirname(current_path))
            status_var.set(f"Saved copy -> {new_fname}")
        return True
    elif choice.get() == "discard":
        return True
    else:   # cancel
        # restore listbox highlight to the current file
        if current_path:
            cur_fname = os.path.basename(current_path)
            items = list(listbox.get(0, tk.END))
            if cur_fname in items:
                idx = items.index(cur_fname)
                listbox.selection_clear(0, tk.END)
                listbox.selection_set(idx)
                listbox.see(idx)
        return False


def on_listbox_select(event=None):
    selection = listbox.curselection()
    if not selection:
        return
    filename = listbox.get(selection[0])
    new_path = os.path.join(folder_entry.cget("text"), filename)
    if new_path == current_path:
        return   # clicked the already-open file
    if not _confirm_switch(new_path):
        return   # user cancelled
    open_file(new_path)
 
 
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
        save_as_btn.config(state=tk.NORMAL)
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
    global current_path
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
    current_path = path
    clear_edited()
    load_files(os.path.dirname(path))
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
 
 
# ── copy / paste helpers ──────────────────────────────────────────────────────

def _top_level_key_of(item):
    """
    Return the top-level JSON key name for any tree item.
    E.g. clicking a leaf deep inside TreeClasses returns 'TreeClasses'.
    Returns None if the item is not rooted at a top-level key.
    """
    current = item
    while current:
        parent = tree.parent(current)
        if parent == "":          # direct child of invisible root → top-level key
            text = tree.item(current, "text")
            # strip branch suffix ("  {}" / "  [N]") or leaf suffix (":  value")
            if ":  " in text:
                return text.split(":  ")[0].strip()
            elif "  " in text:
                return text.split("  ")[0].strip()
            return text.strip()
        current = parent
    return None


def copy_item():
    """
    Copy the right-clicked item to the internal clipboard.

    Two cases:
      • Item is inside TreeClasses → copy that single tree-class dict.
      • Otherwise → copy the entire top-level key/value pair.

    Updates the paste button label to reflect what is in the buffer.
    """
    import copy
    global _clipboard
    item = tree_menu._target_item
    if item is None or current_data is None:
        return

    top_key = _top_level_key_of(item)
    if top_key is None:
        return

    child, tc_idx = _tc_index(item)

    if tc_idx is not None:
        # inside TreeClasses: copy the specific tree-class dict
        tc = current_data["TreeClasses"][tc_idx]
        desc = tc.get("Description", f"class {tc_idx}")
        _clipboard = {
            "key":           "TreeClasses",
            "value":         copy.deepcopy(tc),
            "is_tree_class": True,
            "label":         f"tree class '{desc}'",
        }
    else:
        # any other top-level key
        val = current_data.get(top_key)
        _clipboard = {
            "key":           top_key,
            "value":         copy.deepcopy(val),
            "is_tree_class": False,
            "label":         top_key,
        }

    _update_paste_button()
    status_var.set(f"Copied: {_clipboard['label']}")


def _update_paste_button():
    """Relabel and enable/disable the paste button based on clipboard state."""
    if _clipboard is None:
        paste_btn.config(state=tk.DISABLED, text="Paste")
    else:
        paste_btn.config(state=tk.NORMAL,
                         text=f"Paste {_clipboard['label']}")


def on_paste():
    """
    Paste clipboard contents into the currently loaded JSON.

    • is_tree_class=True  → append the tree-class dict to TreeClasses.
    • is_tree_class=False → overwrite (or add) the top-level key in current_data.

    Marks the document as edited and refreshes the tree.
    """
    import copy
    global _clipboard
    if _clipboard is None or current_data is None:
        return

    val   = copy.deepcopy(_clipboard["value"])   # paste a fresh copy each time
    key   = _clipboard["key"]
    label = _clipboard["label"]

    if _clipboard["is_tree_class"]:
        if "TreeClasses" not in current_data or not isinstance(current_data["TreeClasses"], list):
            current_data["TreeClasses"] = []
        current_data["TreeClasses"].append(val)
        new_idx = len(current_data["TreeClasses"]) - 1
        status_var.set(f"Pasted {label} as class {new_idx}.")
    else:
        new_idx = None
        current_data[key] = val
        status_var.set(f"Pasted {label}.")

    load_tree(current_data)
    mark_edited()

    # after pasting a tree class, open TreeClasses and show the new entry
    if new_idx is not None:
        _expand_and_show_description(new_idx)


# ── file list right-click menu ────────────────────────────────────────────────

def _listbox_item_at(event):
    """Return the filename under the cursor, or None."""
    idx = listbox.nearest(event.y)
    if idx < 0 or idx >= listbox.size():
        return None
    # nearest() can return a row even when the click is below all items;
    # bbox() returns empty string when the index is out of the visible area.
    bbox = listbox.bbox(idx)
    if not bbox:
        return None
    # reject clicks below the last item's bottom edge
    _, y0, _, h = bbox
    if event.y > y0 + h:
        return None
    return listbox.get(idx)


def _duplicate_file():
    """Duplicate the right-clicked JSON file with a _copy suffix."""
    fname = file_menu._target_fname
    if not fname:
        return
    folder = folder_entry.cget("text")
    src    = os.path.join(folder, fname)
    stem, ext = os.path.splitext(fname)

    # find a non-colliding name: stem_copy.json, stem_copy2.json, …
    candidate = os.path.join(folder, stem + "_copy" + ext)
    counter   = 2
    while os.path.exists(candidate):
        candidate = os.path.join(folder, f"{stem}_copy{counter}{ext}")
        counter  += 1

    try:
        import shutil
        shutil.copy2(src, candidate)
    except OSError as e:
        status_var.set(f"Error duplicating: {e}")
        return

    load_files(folder)
    new_fname = os.path.basename(candidate)
    # highlight the new file in the listbox
    items = list(listbox.get(0, tk.END))
    if new_fname in items:
        idx = items.index(new_fname)
        listbox.selection_clear(0, tk.END)
        listbox.selection_set(idx)
        listbox.see(idx)
    status_var.set(f"Duplicated → {new_fname}")


def _write_copy(src_path, data):
    """
    Write data as a new _modified[N].json in the same folder as src_path.
    Strips any existing _modifiedN suffix to avoid chaining.
    Returns the new filename on success, None on failure.
    """
    import re as _re
    folder = os.path.dirname(src_path)
    stem, ext = os.path.splitext(os.path.basename(src_path))
    base_stem = _re.sub(r'_modified\d*$', '', stem)

    candidate = os.path.join(folder, base_stem + "_modified" + ext)
    counter = 2
    while os.path.exists(candidate):
        candidate = os.path.join(folder, f"{base_stem}_modified{counter}{ext}")
        counter += 1

    try:
        with open(candidate, "w") as f:
            json.dump(data, f, indent=2)
    except OSError as e:
        status_var.set(f"Error saving copy: {e}")
        return None
    return os.path.basename(candidate)


def _save_copy_file():
    """Called from the file right-click menu."""
    fname = file_menu._target_fname
    if not fname:
        return
    folder = folder_entry.cget("text")
    src = os.path.join(folder, fname)
    if current_path == src and current_data is not None:
        data_to_write = current_data
    else:
        try:
            with open(src) as f:
                data_to_write = json.load(f)
        except OSError as e:
            status_var.set(f"Error reading source: {e}")
            return
    new_fname = _write_copy(src, data_to_write)
    if new_fname:
        load_files(folder)
        items = list(listbox.get(0, tk.END))
        if new_fname in items:
            ni = items.index(new_fname)
            listbox.selection_clear(0, tk.END)
            listbox.selection_set(ni)
            listbox.see(ni)
        status_var.set(f"Saved copy -> {new_fname}")

def _rename_file():
    """
    Inline-rename the right-clicked JSON file.
    Places a tk.Entry widget directly over the listbox row so the user
    can type the new name and confirm with Enter (or cancel with Escape).
    If the file is currently open, current_path is updated to match.
    """
    fname = file_menu._target_fname
    if not fname:
        return
    folder = folder_entry.cget("text")

    # find the listbox index of the target file
    items = list(listbox.get(0, tk.END))
    if fname not in items:
        return
    lb_idx = items.index(fname)

    # get pixel bbox of that row inside the listbox widget
    bbox = listbox.bbox(lb_idx)
    if not bbox:
        return
    x, y, w, h = bbox

    stem, ext = os.path.splitext(fname)
    entry = tk.Entry(listbox, font=listbox.cget("font") or ("TkDefaultFont", 10))
    entry.place(x=0, y=y, width=listbox.winfo_width(), height=h)
    entry.insert(0, stem)        # pre-fill with current stem (no extension)
    entry.select_range(0, tk.END)
    entry.focus_set()

    def _commit_rename(event=None):
        global current_path
        new_stem = entry.get().strip()
        entry.destroy()
        if not new_stem or new_stem == stem:
            return   # empty or unchanged → cancel
        new_fname = new_stem + ext
        src  = os.path.join(folder, fname)
        dest = os.path.join(folder, new_fname)
        if os.path.exists(dest):
            status_var.set(f"Rename failed: '{new_fname}' already exists.")
            return
        try:
            os.rename(src, dest)
        except OSError as e:
            status_var.set(f"Rename error: {e}")
            return
        # update current_path if the renamed file was open
        if current_path == src:
            current_path = dest
        load_files(folder)
        # re-select renamed file
        new_items = list(listbox.get(0, tk.END))
        if new_fname in new_items:
            ni = new_items.index(new_fname)
            listbox.selection_clear(0, tk.END)
            listbox.selection_set(ni)
            listbox.see(ni)
        status_var.set(f"Renamed → {new_fname}")

    def _cancel_rename(event=None):
        entry.destroy()

    entry.bind("<Return>",  _commit_rename)
    entry.bind("<Escape>",  _cancel_rename)
    entry.bind("<FocusOut>", _commit_rename)


def on_listbox_right_click(event):
    fname = _listbox_item_at(event)
    if not fname:
        return
    # highlight the row under the cursor
    idx = listbox.nearest(event.y)
    listbox.selection_clear(0, tk.END)
    listbox.selection_set(idx)
    file_menu._target_fname = fname
    file_menu.tk_popup(event.x_root, event.y_root)


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
 
# main area — horizontal PanedWindow so the user can drag the file-list width
main = ttk.Frame(root, padding=10)
main.pack(fill=tk.BOTH, expand=True)

h_paned = ttk.PanedWindow(main, orient=tk.HORIZONTAL)
h_paned.pack(fill=tk.BOTH, expand=True)

# left: file list pane
left = ttk.Frame(h_paned)
h_paned.add(left, weight=0)   # weight=0: doesn't steal extra space on resize

ttk.Label(left, text="Files").pack(anchor=tk.W)

lb_frame = ttk.Frame(left)
lb_frame.pack(fill=tk.BOTH, expand=True)

lb_yscroll = ttk.Scrollbar(lb_frame, orient=tk.VERTICAL)
lb_xscroll = ttk.Scrollbar(lb_frame, orient=tk.HORIZONTAL)
listbox = tk.Listbox(lb_frame, activestyle="dotbox",
                     yscrollcommand=lb_yscroll.set,
                     xscrollcommand=lb_xscroll.set)
lb_yscroll.config(command=listbox.yview)
lb_xscroll.config(command=listbox.xview)
lb_xscroll.pack(side=tk.BOTTOM, fill=tk.X)
lb_yscroll.pack(side=tk.RIGHT,  fill=tk.Y)
listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
listbox.bind("<<ListboxSelect>>", on_listbox_select)

file_menu = tk.Menu(root, tearoff=0)
file_menu._target_fname = None
file_menu.add_command(label="Duplicate file",  command=_duplicate_file)
file_menu.add_command(label="Rename file",     command=_rename_file)
listbox.bind("<Button-3>", on_listbox_right_click)

# right column
right = ttk.Frame(h_paned)
h_paned.add(right, weight=1)  # weight=1: right pane gets all extra space
 
# button bar — packed first so it is never squeezed out
btn_bar = ttk.Frame(right)
btn_bar.pack(side=tk.BOTTOM, fill=tk.X, pady=6)
run_btn = ttk.Button(btn_bar, text="▶  Run FRT", command=on_run, state=tk.DISABLED)
run_btn.pack(side=tk.LEFT, padx=(0, 6))
save_as_btn = ttk.Button(btn_bar, text="Save as…", command=on_save_as, state=tk.DISABLED)
save_as_btn.pack(side=tk.LEFT, padx=(0, 6))
paste_btn = ttk.Button(btn_bar, text="Paste", command=on_paste, state=tk.DISABLED)
paste_btn.pack(side=tk.LEFT, padx=(0, 6))
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
    # re-open the TreeClasses branch so it doesn't collapse on the user
    for top_item in tree.get_children(""):
        if tree.item(top_item, "text").startswith("TreeClasses"):
            tree.item(top_item, open=True)
            # scroll to the branch so it stays visible
            tree.see(top_item)
            break
 
 
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
 
 
def _expand_and_show_description(tc_idx):
    """
    After load_tree(), find the tree-class node at tc_idx inside TreeClasses,
    expand it, locate the Description leaf, and scroll it into view.
    """
    # find the TreeClasses branch at the root level
    for top_item in tree.get_children(""):
        text = tree.item(top_item, "text")
        if text.startswith("TreeClasses"):
            # its direct children are the individual tree-class nodes ("0 {}", "1 {}", …)
            tc_children = tree.get_children(top_item)
            if tc_idx < len(tc_children):
                tc_node = tc_children[tc_idx]
                # expand the TreeClasses branch and the specific class node
                tree.item(top_item, open=True)
                tree.item(tc_node, open=True)
                # find the Description leaf among the class's children
                for leaf in tree.get_children(tc_node):
                    if tree.item(leaf, "text").startswith("Description"):
                        tree.see(leaf)
                        tree.selection_set(leaf)
                        return
                # fallback: just scroll the class node into view
                tree.see(tc_node)
            return


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
 
    # expand the node so it is visible and bbox() returns valid coordinates
    tc_parent = tree.parent(child)
    if tc_parent:
        tree.item(tc_parent, open=True)
    tree.see(child)

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
        # expand the tree class node and scroll the Description leaf into view
        _expand_and_show_description(tc_idx)
 
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
tree_menu.add_command(label="Copy",                      command=copy_item)
tree_menu.add_separator()
tree_menu.add_command(label="Delete this tree class",    command=delete_tree_class)
tree_menu.add_command(label="Duplicate this tree class", command=duplicate_tree_class)
tree_menu.add_command(label="Add Description",           command=add_description)
tree_menu.add_separator()
tree_menu.add_command(label="Replace list with file name…",
                      command=lambda: replace_list_with_file(tree_menu._target_item))


def on_tree_right_click(event):
    """
    Called when the user right-clicks anywhere on the tree.

    Shows a context menu for any tree item that has a top-level key:
      • "Copy"                     → always available for any top-level item.
      • Tree-class actions         → only inside TreeClasses.
      • "Replace list with file…"  → only for list [N] branch nodes.
    """
    item = tree.identify_row(event.y)
    if not item:
        return

    tree.selection_set(item)
    tree_menu._target_item = item

    top_key  = _top_level_key_of(item)
    child, idx = _tc_index(item)
    is_tc    = idx is not None
    is_list  = is_list_node(item) is not None
    has_desc = is_tc and _tc_has_description(idx)

    # "Copy" is available whenever the item belongs to a top-level key
    has_top = top_key is not None and current_data is not None
    tree_menu.entryconfig("Copy",
                          state=tk.NORMAL if has_top else tk.DISABLED)

    # tree-class actions
    tree_menu.entryconfig("Delete this tree class",
                          state=tk.NORMAL if is_tc else tk.DISABLED)
    tree_menu.entryconfig("Duplicate this tree class",
                          state=tk.NORMAL if is_tc else tk.DISABLED)
    tree_menu.entryconfig("Add Description",
                          state=tk.NORMAL if (is_tc and not has_desc) else tk.DISABLED)

    # list replacement
    tree_menu.entryconfig("Replace list with file name…",
                          state=tk.NORMAL if is_list else tk.DISABLED)

    # show menu only if at least Copy or one other option is applicable
    if has_top or is_tc or is_list:
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
    total   = paned.winfo_height()
    total_w = h_paned.winfo_width()
    if total > 10 and total_w > 10:
        paned.sashpos(0, int(total * 0.60))
        h_paned.sashpos(0, 220)   # initial file-list width in pixels
    else:
        root.after(100, _set_sash)   # retry if not drawn yet
 
root.after(200, _set_sash)
root.mainloop()
