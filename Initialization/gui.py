#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 16:22:54 2026

@author: oliverlemke
"""

#### Do not run within spyder!!!

### TODO
###### For normalization on peptide list add a dropdown and put in Custom
###### Directory with predefined lists (e.g. OSPP) -> User can add custom.txt
###### Maybe generate a file with all paths to be copied by NextFlow and adjust path if Docker is checked.

##### Layout
#----------------------------------#
#-   flag_1    - flag_parameter_1 -#
#              - flag_parameter_2 -#
#              - flag_parameter_3 -#
#              - flag_parameter_4 -#
#----------------------------------#
#-   flag_2    - flag_parameter_1 -#
#              - flag_parameter_2 -#
#              - flag_parameter_3 -#
#----------------------------------#

import sys
import os
import json
import numpy as np

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget,
    QVBoxLayout, QHBoxLayout, QGroupBox,
    QListWidget, QPushButton, QLabel,
    QLineEdit, QFileDialog, QListWidgetItem,
    QAbstractItemView, QGridLayout, QProgressBar,
    QCheckBox, QComboBox, QFrame, QSizePolicy,
    QStackedWidget, QScrollArea
)
from PySide6.QtCore import Qt, Signal, QTimer
from PySide6.QtGui import QIcon, QIntValidator, QDoubleValidator
from shiboken6 import isValid

import gui_funcs

# ----------------------------
# Data Model
# ----------------------------
class WorkflowStep:
    def __init__(self, name, params=None):
        self.name = name
        # persistent parameter storage
        self.params = {}

# ----------------------------
# Parameter Widgets
# ----------------------------
class DynamicParameterConfig(QWidget):
    id_columns = ["Precursor.Id","Stripped.Sequence","Protein.Group","Genes","Peptide.Id","Proteotypic"]
    def __init__(self, step, flags_dict, params_dict,
                 meta_columns=None,
                 data_columns=None):
        super().__init__()
        self.step = step
        self.flags_dict = {k.upper(): v for k, v in flags_dict.items()}
        self.params_dict = params_dict
        self.meta_columns = meta_columns or []
        self.meta_columns_full = meta_columns if isinstance(meta_columns, dict) else {}
        self.data_columns = data_columns or []
        self.param_widgets = {}
        self.param_widgets_by_name = {}
        self._sync_guard = False
        self.visibility_rules = {}
        self.visibility_targets = {}
        self.meta_lookup = {}

        self.main_layout = QVBoxLayout(self)
        
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Maximum)
        
        self.build_ui()

    # ----------------------------
    # MAIN UI
    # ----------------------------
    def build_ui(self):
        step_flags = self.flags_dict.get(self.step.name.upper(), {})
    
        # Special layout handling
        is_general = self.step.name == "General"
        processed_children = set()

        priority = []
        normal = []
        
        # Build flag page
        for flag_name, meta in step_flags.items():
        
            if isinstance(meta, dict) and meta.get("PriorityFirst", False):
                priority.append((flag_name, meta))
            else:
                normal.append((flag_name, meta))
        
        for flag_name, meta in priority + normal:
            if flag_name in processed_children:
                continue
            box = QGroupBox(flag_name)
            grid = QGridLayout()
            
            (row_widget, desc, flag_widget), child_name = self.build_row_from_meta(
                flag_name,
                meta,
                step_flags
            )
            
            # skip children already rendered
            if child_name:
                processed_children.update(child_name)
            
            # remove left-side label for flags
            layout = row_widget.layout()
            
            if layout and layout.count() > 0:
                first = layout.itemAt(0).widget()
            
                if isinstance(first, QLabel):
                    first.hide()
            
            left_widget = row_widget

            # Parameters container
            param_container = AutoResizingStackedWidget()

            # Build pages for each flag value
            # Empty page first
            has_empty = isinstance(flag_widget, QComboBox) and flag_widget.findText("") != -1
            if has_empty:
                empty_page = QWidget()
                param_container.addWidget(empty_page)
            
            flag_values = self.get_flag_values(flag_widget, meta)
            
            for value in flag_values:
                page = self.build_param_page(flag_name, value)
                param_container.addWidget(page)
            
            # Set initial page
            if isinstance(flag_widget, QCheckBox):
                param_container.setCurrentIndex(1 if flag_widget.isChecked() else 0)
            
            elif isinstance(flag_widget, QComboBox):

                current_index = flag_widget.currentIndex()
            
                if has_empty:
                    # page 0 = empty page
                    param_container.setCurrentIndex(current_index + 1)
                else:
                    param_container.setCurrentIndex(current_index)
            
            # connect change
            if isinstance(flag_widget, QComboBox):

                def update_page(index, pc=param_container):
                    if index < 0:
                        pc.setCurrentIndex(0)
                    elif not has_empty:
                        pc.setCurrentIndex(index)
                    else:
                        pc.setCurrentIndex(index + 1)
            
                flag_widget.currentIndexChanged.connect(update_page)
            
            elif isinstance(flag_widget, QCheckBox):
            
                flag_widget.toggled.connect(
                    lambda checked, pc=param_container:
                        pc.setCurrentIndex(1 if checked else 0)
                )
                        
            param_container.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Maximum)
            
            # Define flag - parameter separator
            separator = QFrame()
            separator.setFrameShape(QFrame.VLine)
            separator.setFrameShadow(QFrame.Plain)
            separator.setLineWidth(1)
            separator.setStyleSheet("""
                QFrame {
                    color: #e6e6e6;
                    margin-left: 6px;
                    margin-right: 6px;
                }
            """)
            separator.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)
            
            # Layout proportions
            if is_general:
                grid.setColumnStretch(0, 9)
                grid.setColumnStretch(1, 0)
                grid.setColumnStretch(2, 1)
            else:
                grid.setColumnStretch(0, 1)
                grid.setColumnStretch(1, 0)
                grid.setColumnStretch(2, 3)

            grid.addWidget(left_widget, 0, 0)
            grid.addWidget(desc, 1, 0)
            
            # For General: no right-side parameter area
            if not is_general:
                grid.addWidget(separator, 0, 1, 2, 1)
                grid.addWidget(param_container, 0, 2, 2, 1)

            box.setLayout(grid)
            box.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Maximum)
            self.main_layout.addWidget(box)
            
        if not hasattr(self, "_visibility_connected"):
            self._visibility_connected = True
        
            for w in self.param_widgets.values():
                if isinstance(w, QComboBox):
                    w.currentIndexChanged.connect(lambda *_: self.update_visibility())
                elif isinstance(w, QCheckBox):
                    w.toggled.connect(lambda *_: self.update_visibility())
        self.update_visibility()

    # ----------------------------
    # ROW CONSTRUCTOR
    # ----------------------------
    def build_row_from_meta(self, name, meta, source_dict):
    
        child_names = meta.get("Child")

        if isinstance(child_names, str):
            child_names = [child_names]
    
        # Dependency case
        if child_names:

            child_metas = [
                source_dict[c]
                for c in child_names
                if c in source_dict
            ]
        
            if child_metas:
        
                row, desc = self.build_dependency_row(
                    name,
                    meta,
                    child_names,
                    child_metas
                )
        
                return (row, desc, self.param_widgets[name]), child_names
    
        # Standard row
        row, desc = self.build_param_row(name, meta)

        return (row, desc, self.param_widgets[name]), None

    # ----------------------------
    # PARAM PAGE
    # ----------------------------
    def build_param_page(self, flag_name, flag_value):
        page = QWidget()
        page.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Maximum)
        layout = QVBoxLayout(page)
    
        params = self.params_dict.get(flag_name, {}).get(flag_value, {})
    
        # Safety check
        if not isinstance(params, dict):
            return page
    
        processed_children = set()
    
        for param_name, meta in params.items():
    
            # Safety check
            if not isinstance(meta, dict):
                continue
    
            if param_name in processed_children:
                continue
    
            (row, desc, _), child_name = self.build_row_from_meta(
                param_name,
                meta,
                params
            )
    
            if child_name:
                processed_children.update(child_name)
    
            layout.addWidget(row)
            layout.addWidget(desc)
    
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(2)
        layout.setAlignment(Qt.AlignTop)
    
        return page
    
    
    # ----------------------------
    # PARAM ROW 
    # ----------------------------
    def build_param_row(self, name, meta):
        row = QWidget()
        layout = QHBoxLayout(row)

        label = QLabel(name)
        self.style_required_label(label, meta)
        self.meta_lookup[name] = meta

        parsed_type = self.parse_gui_type(meta["GUI"], meta.get("Options"))
        widget = self.create_widget(parsed_type, meta["GUI"])
        
        self.apply_required_validator(
            widget,
            meta.get("NoneValid", False)
        )

        if meta["GUI"].startswith("R"):
            self.apply_numeric_range_validator(
                widget,
                meta.get("Options")
            )
            
        self.register_visibility_dependency(
            name,
            meta.get("VisibleIf")
        )

        # store reference
        self.param_widgets[name] = widget
        self.visibility_targets[name] = row
        self.param_widgets_by_name.setdefault(name, []).append(widget)
        self.connect_widget_signal(name, widget)
        
        if meta["GUI"].startswith("S2"):
            self._link_min_max(widget)

        # Restore saved value
        if name in self.step.params:
        
            self.apply_default(
                widget,
                self.step.params[name]
            )
        else:
        
            self.apply_default(
                widget,
                meta.get("Default")
            )
            
        self.save_param_value(name, widget)

        info_btn, desc = self.create_info(meta.get("Comment", ""))
        desc.setIndent(24)

        layout.addWidget(label)
        layout.addWidget(widget, 1)
        layout.addWidget(info_btn)
        label.setAttribute(Qt.WA_LayoutUsesWidgetRect, True)
        label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)

        return row, desc
    
    # ----------------------------
    # DEPENDENCY PARAM ROW 
    # ----------------------------
    def build_dependency_row(self,
                         parent_name,
                         parent_meta,
                         child_names, 
                         child_metas):

        row = QWidget()
    
        layout = QHBoxLayout(row)
    
        label = QLabel(parent_name)
        self.style_required_label(label, parent_meta)
        self.meta_lookup[parent_name] = parent_meta
    
        parent_type = self.parse_gui_type(
            parent_meta["GUI"],
            parent_meta.get("Options")
        )
    
        # Multi dependency and empty
        allow_empty = "E" in parent_meta["GUI"]

        if "A" in parent_meta["GUI"]:

            widget_to_add = MultiDependencyWidget(
                metadata_dict=self.meta_columns_full,
                data_list=self.data_columns,
                parent_type=parent_type,
                child_metas=child_metas,
                create_widget_func=self.create_widget,
                parse_gui_type_func=self.parse_gui_type,
                allow_empty=allow_empty
            )
            self.apply_required_validator(
                widget_to_add,
                parent_meta.get("NoneValid", False)
            )
        
            # store controlling widget
            self.param_widgets[parent_name] = widget_to_add
            self.visibility_targets[parent_name] = row
            self.param_widgets_by_name.setdefault(parent_name, []).append(widget_to_add)
        
        # Single dependency
        else:
        
            parent_widget = self.create_widget(
                parent_type,
                parent_meta["GUI"]
            )

            if parent_meta["GUI"].startswith("R"):
                self.apply_numeric_range_validator(
                    parent_widget,
                    parent_meta.get("Options")
                )
        
            self.apply_default(
                parent_widget,
                parent_meta.get("Default")
            )
        
            child_widgets = []

            child_container = QWidget()

            child_layout = QHBoxLayout(child_container)
            child_layout.setContentsMargins(0, 0, 0, 0)
            
            child_layout.addWidget(parent_widget, 1)
            
            for child_meta in child_metas:
            
                child_type = self.parse_gui_type(
                    child_meta["GUI"],
                    child_meta.get("Options")
                )
            
                child_widget = self.create_widget(
                    child_type,
                    child_meta["GUI"]
                )
                
                self.apply_default(
                    child_widget,
                    child_meta.get("Default")
                )
            
                dep_widget = DependencyWidget(
                    parent_widget=parent_widget,
                    child_widget=child_widget,
                    metadata_dict=self.meta_columns_full,
                    child_source=child_meta.get("M_column"),
                    use_dependency="X" in child_meta["GUI"],
                    allow_empty=allow_empty
                )
                self.connect_widget_signal(parent_name, dep_widget)
            
                child_layout.addWidget(child_widget, 1)
            
                child_widgets.append(dep_widget)
                
            child_container.child_widgets = child_widgets
            widget_to_add = child_container
            
            self.apply_required_validator(
                widget_to_add,
                parent_meta.get("NoneValid", False)
            )
        
            self.param_widgets[parent_name] = widget_to_add
            self.visibility_targets[parent_name] = row
            
        self.connect_widget_signal(parent_name, widget_to_add)

        # Restore saved value or apply default
        saved_value = self.step.params.get(parent_name, None)
        
        if saved_value is not None:
            self.apply_default(widget_to_add, saved_value)
        else:
            self.apply_default(widget_to_add, parent_meta.get("Default"))
                    
        self.save_param_value(parent_name, widget_to_add)
        
        # Info button
        info_btn, desc = self.create_info(
            parent_meta.get("Comment", "")
        )
        desc.setIndent(24)
    
        layout.addWidget(label)
        layout.addWidget(widget_to_add, 1)
        layout.addWidget(info_btn)
    
        return row, desc

    # ----------------------------
    # HELPERS
    # ----------------------------
    # Info button
    def create_info(self, text):
        info_button = QPushButton()
    
        icon = QIcon.fromTheme("help-about")
        if icon.isNull():
            info_button.setText("?")
        else:
            info_button.setIcon(icon)
    
        info_button.setFixedSize(20, 20)
        info_button.setFlat(True)
        info_button.setCursor(Qt.PointingHandCursor)
        
        info_button.setStyleSheet("border: none;")
    
        label = QLabel(text)
        label.setWordWrap(True)
        label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        label.setMaximumHeight(60)
        label.setStyleSheet("color: gray; font-size: 11px;")
        label.setVisible(False)
    
        info_button.clicked.connect(
            lambda: label.setVisible(not label.isVisible())
        )
    
        return info_button, label
    
    # Extract flags
    def get_flag_values(self, widget, meta):
        if isinstance(widget, QComboBox):
            return [widget.itemText(i) for i in range(widget.count())]

        elif isinstance(widget, QCheckBox):
            return ["0", "1"]

        return ["default"]
    
    # Set default or restore values
    def apply_default(self, widget, value):
    
        if isinstance(value, str):
        
            # tuple syntax: "(0,2)"
            if value.startswith("("):
        
                try:
                    v = value.strip("()").split(",")
        
                    parsed = []
        
                    for x in v:
        
                        x = x.strip()
        
                        if not x:
                            parsed.append(None)
        
                        elif "inf" in x:
                            parsed.append(np.inf)
        
                        elif "." in x:
                            parsed.append(float(x))
        
                        else:
                            parsed.append(int(x))
        
                    value = tuple(parsed)
        
                except:
                    pass
        
            # list syntax: "[1,-5,3]"
            elif value.startswith("["):
        
                try:
                    v = value.strip("[]").split(",")
        
                    parsed = []
        
                    for x in v:
        
                        x = x.strip()
        
                        if not x:
                            parsed.append(None)
        
                        elif "." in x:
                            parsed.append(float(x))
        
                        else:
                            parsed.append(int(x))
        
                    value = parsed
        
                except:
                    pass        
                
        # Dependency restore
        if hasattr(widget, "child_widgets"):
        
            if value is None:
                return
        
            if not isinstance(value, list):
                value = [value]
        
            for dep, val in zip(widget.child_widgets, value):
                dep.set_value(val)
        
            return
        
        # Multi value restore
        if isinstance(widget, MultiValueWidget):
    
            if value is None or value == "None":
                values = [None]
    
            elif isinstance(value, list):
                values = value
    
            else:
                values = [value]
    
            while len(widget.rows) < len(values):
                widget.add_row()
    
            for i, val in enumerate(values):
    
                row_widget, _, _ = widget.rows[i]
    
                field = row_widget.layout().itemAt(0).widget()
    
                if isinstance(field, BaseParamWidget):
                    field.set_value(val)
    
            while len(widget.rows) > len(values):
                row_widget, _, _ = widget.rows[-1]
                widget.remove_row(row_widget)
    
            return
    
        # Adapter handling
        if hasattr(widget, "set_value") and hasattr(widget, "get_value"):
        
            widget.blockSignals(True)
        
            if value is None or value == "None":
                widget.clear_value()
            else:
                widget.set_value(value)
        
            widget.blockSignals(False)
        
            self.revalidate_widget(widget)
    
    # Check the validity of the values
    def revalidate_widget(self, widget):
    
        if hasattr(widget, "_validate_required"):
            widget._validate_required()
        if hasattr(widget, "textChanged"):
            try:
                widget.textChanged.emit(widget.text())
            except:
                pass
    
        if hasattr(widget, "currentTextChanged"):
            try:
                widget.currentTextChanged.emit(widget.currentText())
            except:
                pass
    
        if hasattr(widget, "value_changed"):
            try:
                widget.value_changed.emit()
            except:
                pass
    
    # link fields in tuple
    def _link_min_max(self, widget):
        if not isinstance(widget, QWidget):
            return
    
        layout = widget.layout()
        if not layout or layout.count() < 2:
            return
    
        min_field = widget.min_field
        max_field = widget.max_field
    
        if not isinstance(min_field, QLineEdit) or not isinstance(max_field, QLineEdit):
            return
    
        def enforce():
            min_text = min_field.text().strip()
            max_text = max_field.text().strip()
        
            if not min_text or not max_text:
                return
        
            try:
                min_val = float(min_text)
                max_val = float(max_text)
        
                if max_val < min_val:
                    max_field.setText(min_field.text())
            except:
                pass
        min_field.textChanged.connect(enforce)
        max_field.textChanged.connect(enforce)
        
        def read_range(self, widget):
            layout = widget.layout()
            if not layout or layout.count() < 2:
                return None
        
            min_field = layout.itemAt(0).widget()
            max_field = layout.itemAt(1).widget()
        
            min_text = min_field.text().strip()
            max_text = max_field.text().strip()
        
            min_val = float(min_text) if min_text else None
            max_val = float(max_text) if max_text else None
        
            if min_val is None and max_val is None:
                return None
        
            if max_val is None:
                max_val = np.inf
        
            return (min_val, max_val)
            
    # ----------------------------
    # GUI TYPE PARSER
    # ----------------------------
    def parse_gui_type(self, gui, options=None):
        if not gui:
            return ("str", False)
    
        multi = "A" in gui
        code = gui
    
        if code.startswith("M"):
            return ("meta", multi)
    
        if code.startswith("D"):
            return ("data", multi)
    
        if code.startswith("R1"):
            return ("float", multi)
    
        if code.startswith("R3"):
            return ("int", multi)
    
        if code.startswith("T"):
            return ("str", multi)
    
        elif code.startswith("R"):
            return ("float", multi)
    
        elif code.startswith("C"):
            return ("bool", False)
    
        elif code.startswith("L"):
            return (options or [], multi)
    
        elif code.startswith("S2"):
            return (("min", "max"), multi)
        
        elif code.startswith("S3"):
            return (("v1", "v2", "v3"), multi)
        
        elif code.startswith("F"):
            return ("file", multi)
        
        elif gui.startswith("I"):
            return ("id_columns", multi)

    # Validator for values that can not be empty
    def apply_required_validator(self, widget, none_valid=False):

        if not hasattr(self, "_validator_guard"):
            self._validator_guard = set()
    
        def validate():
    
            key = id(widget)
    
            if key in self._validator_guard:
                return
    
            self._validator_guard.add(key)
    
            try:
    
                valid = True
    
                # ----------------------------
                # QLineEdit
                # ----------------------------
                if isinstance(widget, QLineEdit):
    
                    valid = bool(widget.text().strip())
                    
                # ----------------------------
                # File Browser
                # ----------------------------                    
                elif isinstance(widget, FileBrowserAdapter):

                    valid = bool(widget.get_value())
    
                # ----------------------------
                # QComboBox
                # ----------------------------
                elif isinstance(widget, QComboBox):
    
                    valid = bool(widget.currentText().strip())
    
                # ----------------------------
                # MultiValueWidget
                # ----------------------------
                elif isinstance(widget, MultiValueWidget):
    
                    for i, (row_widget, _, _) in enumerate(widget.rows):
    
                        field = row_widget.layout().itemAt(0).widget()
    
                        if i == 0:
    
                            if isinstance(field, QLineEdit):
                                valid &= bool(field.text().strip())
    
                            elif isinstance(field, QComboBox):
                                valid &= bool(field.currentText().strip())
    
                            elif hasattr(field, "get_value"):
                                valid &= field.get_value() not in [None, "", [], ()]
    
                # ----------------------------
                # DependencyWidget
                # ----------------------------
                elif isinstance(widget, DependencyWidget):
    
                    parent = widget.parent_widget
                    child = widget.child_widget
    
                    parent_val = parent.currentText().strip()
    
                    if isinstance(child, QComboBox):
                        child_val = child.currentText().strip()
    
                    elif isinstance(child, QLineEdit):
                        child_val = child.text().strip()
    
                    elif hasattr(child, "get_value"):
                        child_val = child.get_value()
    
                    else:
                        child_val = None
    
                    valid = bool(parent_val) and bool(child_val)
    
                # ----------------------------
                # MultiDependencyWidget
                # ----------------------------
                elif isinstance(widget, MultiDependencyWidget):
                    valid = True
                    if widget.rows:
                        row_widget, _, _ = widget.rows[0]
                        layout = row_widget.layout()
                        if layout is not None:
                            dep = layout.itemAt(0).widget()
                            if dep is not None:
                                parent_val = dep.parent_combo.get_value()
                                if not parent_val:
                                    valid = False
                                else:
                                    for child_widget, _ in dep.child_widgets:
                                        if hasattr(child_widget, "get_value"):
                                            child_val = child_widget.get_value()
                                        elif isinstance(child_widget, QComboBox):
                                            child_val = child_widget.currentText()
                                        elif isinstance(child_widget, QLineEdit):
                                            child_val = child_widget.text()
                                        else:
                                            child_val = None
                                        if child_val in [None, ""]:
                                            valid = False
                                            break
    
                # ----------------------------
                # Generic BaseParamWidget
                # ----------------------------
                elif hasattr(widget, "get_value"):
    
                    value = widget.get_value()

                    if hasattr(widget, "is_valid"):
                        valid = widget.is_valid()
                    else:
                        valid = value not in [None, "", [], ()]
    
                # NoneValid override
                if none_valid:
                    return
    
                self.apply_error(widget, not valid)
    
            finally:
    
                self._validator_guard.discard(key)
    
        widget._validate_required = validate
    
        # ----------------------------
        # CONNECT SIGNALS
        # ----------------------------
        if isinstance(widget, QLineEdit):
    
            widget.textChanged.connect(validate)
    
        elif isinstance(widget, QComboBox):
    
            widget.currentTextChanged.connect(validate)
    
        elif isinstance(widget, QCheckBox):
    
            widget.stateChanged.connect(validate)
            
        elif isinstance(widget, FileBrowserAdapter):

            widget.line.textChanged.connect(validate)
    
        elif hasattr(widget, "value_changed"):
    
            widget.value_changed.connect(validate)
    
        validate()

    # Validator for numeric values
    def apply_numeric_range_validator(self, widget, options):

        if not options:
            return
    
        if isinstance(options, list):
            options = options[0]
    
        if not isinstance(options, str):
            return
    
        text = options.strip()
    
        if len(text) < 5:
            return
    
        left_inclusive = text.startswith("[")
        right_inclusive = text.endswith("]")
    
        try:
            body = text[1:-1]
            min_str, max_str = [x.strip() for x in body.split(",")]
    
            min_val = float(min_str) if min_str else -np.inf
            max_val = float(max_str) if max_str else np.inf
    
        except:
            return
    
        def validate():
    
            value_text = widget.text().strip()
    
            if not value_text:
                return
    
            try:
                value = float(value_text)
    
            except:
                return
    
            invalid = False
    
            if left_inclusive:
                invalid |= value < min_val
            else:
                invalid |= value <= min_val
    
            if right_inclusive:
                invalid |= value > max_val
            else:
                invalid |= value >= max_val
    
            self.apply_error(widget, invalid)
    
        widget.textChanged.connect(validate)
        
    # ----------------------------
    # WIDGET FACTORY
    # ----------------------------
    def create_widget(self, var_type, gui=None):
        if isinstance(var_type, tuple) and len(var_type) == 2 and isinstance(var_type[1], bool):
            base_type, multi = var_type
        else:
            base_type = var_type
            multi = False
    
        def make_single():
            return self._create_single_widget(base_type, gui)
    
        if multi:
            return MultiValueWidget(make_single)
        else:
            return make_single()
        
    def _create_single_widget(self, var_type, gui=None):

        allow_empty = gui and "E" in gui
    
        # ----------------------------
        # META
        # ----------------------------
        if var_type == "meta":
    
            combo = ComboBoxAdapter()
    
            if allow_empty:
                combo.addItem("")
    
            keys = [k for k in self.meta_columns_full.keys()
                    if k != "__dummy_type__"]
    
            combo.addItems(keys)
            combo.setCurrentIndex(-1)
    
            return combo
    
        # ----------------------------
        # DATA
        # ----------------------------
        if var_type == "data":
    
            combo = ComboBoxAdapter()
    
            if allow_empty:
                combo.addItem("")
    
            if gui and "D1P" in gui:
                combo.addItem("Peptide.Id")
    
            combo.addItems(self.data_columns)
    
            combo.setCurrentIndex(0 if allow_empty else -1)
    
            return combo
    
        # ----------------------------
        # RANGE
        # ----------------------------
        if isinstance(var_type, tuple):

            # S2
            if var_type == ("min", "max"):
                integer = gui and "3" in gui
                return RangeWidgetAdapter(integer=integer)
        
            # S3
            elif var_type == ("v1", "v2", "v3"):
                widget = TripleIntWidget()   # create this class
                return widget
    
        # ----------------------------
        # BOOL
        # ----------------------------
        if var_type == "bool":
            return CheckBoxAdapter()
    
        # ----------------------------
        # LIST
        # ----------------------------
        elif isinstance(var_type, list):
    
            combo = ComboBoxAdapter()
    
            combo.addItems([str(v) for v in var_type])
            combo.setCurrentIndex(-1)
    
            return combo
    
        # ----------------------------
        # INT
        # ----------------------------
        elif var_type == "int":
    
            line = LineEditAdapter()
            line.setValidator(QIntValidator())
    
            return line
    
        # ----------------------------
        # FLOAT
        # ----------------------------
        elif var_type == "float":
    
            line = LineEditAdapter()
            line.setValidator(QDoubleValidator())
    
            return line
        
        # ----------------------------
        # FILE BROWSER
        # ----------------------------
        elif var_type == "file":
            return FileBrowserAdapter()
        
        # ----------------------------
        # ID COLUMNS
        # ----------------------------
        elif var_type == "id_columns":

            combo = ComboBoxAdapter()
            combo.addItems(DynamicParameterConfig.id_columns)
            combo.setCurrentIndex(-1)
        
            return combo
        
        # ----------------------------
        # DEFAULT STRING
        # ----------------------------
        return LineEditAdapter()
    
    def save_param_value(self, name, widget):
        #print(f"[SAVE] {name} -> {type(widget).__name__}")
        
        if hasattr(widget, "child_widgets"):

            val = []
        
            for dep in widget.child_widgets:
                val.append(dep.get_value())
        
            self.step.params[name] = val
        
        elif hasattr(widget, "get_value"):
            val = widget.get_value()
            #print(f"value: {val}")
            self.step.params[name] = val
            if not getattr(widget, "_syncing", False):
                self.sync_linked_widgets(name, widget)
            #print(f"[STORE] {name} = {val}")
        
        # Update shared ID column list
        if name == "input_ids_columns":
        
            values = val if isinstance(val, list) else [val]
        
            # remove old dynamic entries
            defaults = {"Precursor.Id","Stripped.Sequence", "Protein.Group", "Genes", "Peptide.Id", "Proteotypic"}
        
            DynamicParameterConfig.id_columns = [x for x in DynamicParameterConfig.id_columns if x in defaults]
        
            # add new entries
            for v in values:
                if v and v not in DynamicParameterConfig.id_columns:
                    DynamicParameterConfig.id_columns.append(v)
                                
    def connect_widget_signal(self, name, widget):

        # ----------------------------
        # LineEdit
        # ----------------------------
        if isinstance(widget, QLineEdit):
    
            widget.textChanged.connect(
                lambda *_: self.save_param_value(name, widget)
            )
    
        # ----------------------------
        # ComboBox
        # ----------------------------
        elif isinstance(widget, QComboBox):
    
            widget.currentTextChanged.connect(
                lambda *_: self.save_param_value(name, widget)
            )
    
        # ----------------------------
        # Checkbox
        # ----------------------------
        elif isinstance(widget, QCheckBox):
    
            widget.stateChanged.connect(
                lambda *_: self.save_param_value(name, widget)
            )
    
        # ----------------------------
        # Multi widgets
        # ----------------------------
        elif hasattr(widget, "child_widgets"):
        
            for dep in widget.child_widgets:
                dep.value_changed.connect(
                    lambda *_: self.save_param_value(name, widget)
                )
        
        elif isinstance(widget, MultiValueWidget):
    
            widget.value_changed.connect(
                lambda: self.save_param_value(name, widget)
            )
    
        # ----------------------------
        # Range widgets
        # ----------------------------
        elif isinstance(widget, RangeWidgetAdapter):

            def emit():
                self.save_param_value(name, widget)
                self.apply_required_validator(widget)
        
            widget.min_field.textChanged.connect(lambda *_: emit())
            widget.max_field.textChanged.connect(lambda *_: emit())
            
        elif isinstance(widget, DependencyWidget):
            widget.value_changed.connect(
                lambda: self.save_param_value(name, widget)
            )
        
        elif isinstance(widget, DependencyRowWidget):
            widget.value_changed.connect(
                lambda: self.save_param_value(name, widget)
            )
            
        elif isinstance(widget, MultiDependencyWidget):
            widget.value_changed.connect(
                lambda: self.save_param_value(name, widget)
            )
            
        # ----------------------------
        # TripleInt widgets
        # ----------------------------    
        elif isinstance(widget, TripleIntWidget):

            def emit():
                self.save_param_value(name, widget)
                self.apply_required_validator(widget)
        
            widget.f1.textChanged.connect(lambda *_: emit())
            widget.f2.textChanged.connect(lambda *_: emit())
            widget.f3.textChanged.connect(lambda *_: emit())
        
        # ----------------------------
        # File Browser
        # ----------------------------   
        elif isinstance(widget, FileBrowserAdapter):
            widget.line.textChanged.connect(lambda *_: self.save_param_value(name, widget))
            
    def style_required_label(self, label, meta):
    
        none_valid = meta.get("NoneValid", False)
    
        required = none_valid
    
        if required:
    
            label.setStyleSheet("""
                QLabel {
                    color: #999999;
                }
            """)
    
        else:
    
            label.setStyleSheet("")
            
    def sync_linked_widgets(self, name, source_widget):

        if self._sync_guard:
            return
    
        if not hasattr(source_widget, "get_value"):
            return
    
        value = source_widget.get_value()
    
        self._sync_guard = True
    
        try:
    
            for widget in self.param_widgets_by_name.get(name, []):
    
                if widget is source_widget:
                    continue
    
                if not hasattr(widget, "set_value"):
                    continue
    
                # IMPORTANT
                widget._syncing = True
    
                try:
                    widget.set_value(value)
    
                finally:
                    widget._syncing = False
    
            self.step.params[name] = value
    
        finally:
            self._sync_guard = False
            
    def apply_error(self, widget, error):

        targets = self._unwrap_all(widget)
    
        for t in targets:
    
            if not isValid(t):
                continue
    
            try:
                if error:
                    t.setStyleSheet("border: 1px solid red;")
                else:
                    t.setStyleSheet("")
            except RuntimeError:
                pass
    
    def _unwrap(self, widget):

        if isinstance(widget, RangeWidgetAdapter):
            return widget.min_field
    
        if isinstance(widget, DependencyWidget):
            # treat parent as representative (same as non-dependency widgets)
            return widget.parent_widget
    
        if isinstance(widget, DependencyRowWidget):
            return widget.parent_combo
    
        if isinstance(widget, MultiDependencyWidget):
            if widget.rows:
                row_widget = widget.rows[0][0]
                dep = row_widget.layout().itemAt(0).widget()
                return dep.parent_combo  # or dep.child_combo depending on meta rule
    
        if isinstance(widget, MultiValueWidget):
            return widget.rows[0][0].layout().itemAt(0).widget()
    
        return widget
    
    def _unwrap_all(self, widget):

        if isinstance(widget, DependencyWidget):
            return [widget.parent_widget, widget.child_widget]
    
        if isinstance(widget, DependencyRowWidget):
    
            result = [widget.parent_combo]
    
            child = widget.child_widget
    
            if hasattr(child, "get_value"):
                result.append(child)
            elif hasattr(child, "currentText"):
                result.append(child)
            elif hasattr(child, "text"):
                result.append(child)
    
            return result
    
        if isinstance(widget, MultiDependencyWidget):

            if widget.rows:
        
                dep = widget.rows[0][0].layout().itemAt(0).widget()
        
                result = [dep.parent_combo]
        
                for child_widget, _ in dep.child_widgets:
                    result.append(child_widget)
        
                return result
    
        if isinstance(widget, RangeWidgetAdapter):
            return [widget.min_field, widget.max_field]
    
        return [self._unwrap(widget)]
    
    def register_visibility_dependency(self, param_name, rule):
        if not rule:
            return
    
        self.visibility_rules[param_name] = rule
        
    def update_visibility(self):
        for param_name, rule in self.visibility_rules.items():
    
            widget = self.visibility_targets.get(param_name)
            if not widget:
                continue

            visible = True
    
            for controller_name, expected_value in rule.items():
    
                controller = self.param_widgets.get(controller_name)
                if not controller:
                    continue
    
                actual = self.get_widget_value(controller)
    
                if isinstance(expected_value, list):
                    if actual not in expected_value:
                        visible = False
                        break
                else:
                    if actual != expected_value:
                        visible = False
                        break
    
            # apply visibility            
            if param_name in self.flags_dict[self.step.name.upper()]:
                widget.parentWidget().setVisible(visible)
            else:
                widget.setVisible(visible)
                
        self.apply_filtered_options()
    
    def get_widget_value(self, widget):
    
        if isinstance(widget, QComboBox):
            return widget.currentText().strip()
    
        if isinstance(widget, QCheckBox):
            return widget.isChecked()
        
        if hasattr(widget, "get_value"):
            return widget.get_value()
    
        return None
        
    def apply_filtered_options(self):
        for name, widget in self.param_widgets.items():
    
            if not isinstance(widget, QComboBox):
                continue
    
            meta = self.meta_lookup.get(name)
            if not meta:
                continue
    
            rule = meta.get("FilteredOptionsIf")
            if not rule:
                continue
    
            current_value = widget.currentText().strip()
    
            original = meta.get("Options", [])
    
            widget.blockSignals(True)
    
            widget.clear()
            widget.addItems(original)
    
            for controller_name, cases in rule.items():
    
                controller = self.param_widgets.get(controller_name)
                if not controller:
                    continue
    
                actual = self.get_widget_value(controller)
    
                to_remove = cases.get(actual, [])
    
                for item in to_remove:
                    index = widget.findText(item)
                    if index != -1:
                        widget.removeItem(index)
    
            idx = widget.findText(current_value)
            if idx != -1:
                widget.setCurrentIndex(idx)
    
            widget.blockSignals(False)

# ----------------------------
# Automated sizing
# ----------------------------

class AutoResizingStackedWidget(QStackedWidget):
    def sizeHint(self):
        if self.currentWidget():
            return self.currentWidget().sizeHint()
        return super().sizeHint()

    def minimumSizeHint(self):
        if self.currentWidget():
            return self.currentWidget().minimumSizeHint()
        return super().minimumSizeHint()

    def hasHeightForWidth(self):
        return False

    def setCurrentIndex(self, index):
        super().setCurrentIndex(index)
        self.updateGeometry()

# ----------------------------
# Base Adapter Interface
# ----------------------------

class BaseParamWidget:
    """
    Common interface for all parameter widgets.
    """

    def clear_value(self):
        raise NotImplementedError(f"{self.__class__.__name__}.set_value() not implemented")

    def set_value(self, value):
        raise NotImplementedError(f"{self.__class__.__name__}.set_value() not implemented")

    def get_value(self):
        raise NotImplementedError(f"{self.__class__.__name__}.set_value() not implemented")

# ----------------------------
# Multi value Widget
# ----------------------------

class MultiValueWidget(QWidget, BaseParamWidget):
    value_changed = Signal()
    def __init__(self, factory_func):
        super().__init__()
        self.factory_func = factory_func
        self._internal_update = False

        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.setSpacing(2)

        self.rows = []

        self.add_row()

    def add_row(self):

        row_widget = QWidget()
        row_layout = QHBoxLayout(row_widget)
        row_layout.setContentsMargins(0, 0, 0, 0)

        field = self.factory_func()

        if isinstance(field, QComboBox):
            field.currentTextChanged.connect(
                lambda *_: self.value_changed.emit()
            )
        
        elif isinstance(field, QLineEdit):
            field.textChanged.connect(
                lambda *_: self.value_changed.emit()
            )
        
        elif isinstance(field, QCheckBox):
            field.stateChanged.connect(
                lambda *_: self.value_changed.emit()
            )

        add_btn = QPushButton("+")
        add_btn.setFixedWidth(25)

        remove_btn = QPushButton("-")
        remove_btn.setFixedWidth(25)

        # --- buttons container ---
        btn_container = QWidget()
        btn_layout = QHBoxLayout(btn_container)
        btn_layout.setContentsMargins(0, 0, 0, 0)
        btn_layout.setSpacing(2)
        
        btn_layout.addWidget(add_btn)
        btn_layout.addWidget(remove_btn)
        
        btn_container.setFixedWidth(60)  # keeps layout stable
        
        # --- add to row ---
        row_layout.addWidget(field, 1)   # stretch
        row_layout.addWidget(btn_container)

        add_btn.clicked.connect(self.add_row)
        remove_btn.clicked.connect(lambda: self.remove_row(row_widget))

        self.layout.addWidget(row_widget)
        self.rows.append((row_widget, add_btn, remove_btn))

        self.update_buttons()

        self.value_changed.emit()

    def remove_row(self, row_widget):

        if len(self.rows) <= 1:
            return
    
        for i, (row, _, _) in enumerate(self.rows):
    
            if row == row_widget:
    
                self.rows.pop(i)
    
                row.setParent(None)
                row.deleteLater()
    
                break
    
        self.update_buttons()
    
        if not self._internal_update:
            self.value_changed.emit()

    def update_buttons(self):
        n = len(self.rows)

        for i, (row, add_btn, remove_btn) in enumerate(self.rows):

            add_btn.setEnabled(i == n - 1)
            add_btn.setVisible(True)

            remove_btn.setVisible(n > 1)
            
    def is_valid(self):

        for row_widget, _, _ in self.rows:
    
            field = row_widget.layout().itemAt(0).widget()
    
            # BaseParamWidget
            if hasattr(field, "get_value"):
    
                value = field.get_value()
    
                if value in [None, "", [], ()]:
                    return False
    
            # QLineEdit fallback
            elif isinstance(field, QLineEdit):
    
                if not field.text().strip():
                    return False
    
            # QComboBox fallback
            elif isinstance(field, QComboBox):
    
                if not field.currentText().strip():
                    return False
    
        return True
            
    def get_value(self):
    
        values = []
    
        for row_widget, _, _ in self.rows:
    
            field = row_widget.layout().itemAt(0).widget()
    
            if isinstance(field, BaseParamWidget):
    
                values.append(field.get_value())
    
        return values
    
    def set_value(self, values):
    
        if self._internal_update:
            return
    
        self._internal_update = True
    
        try:
    
            if values is None:
                values = [None]
    
            elif not isinstance(values, list):
                values = [values]
    
            # Add missing row
            while len(self.rows) < len(values):
                self.add_row()
    
            # set values
            for i, value in enumerate(values):
    
                row_widget, _, _ = self.rows[i]
    
                field = row_widget.layout().itemAt(0).widget()
    
                # ----------------------------
                # BaseParamWidget
                # ----------------------------
                if isinstance(field, BaseParamWidget):
    
                    field.blockSignals(True)
    
                    try:
    
                        if value in [None, "None"]:
                            field.clear_value()
    
                        else:
                            field.set_value(value)
    
                    finally:
    
                        field.blockSignals(False)
    
                # ----------------------------
                # QLineEdit
                # ----------------------------
                elif isinstance(field, QLineEdit):
    
                    field.blockSignals(True)
    
                    try:
                        field.setText(
                            "" if value is None else str(value)
                        )
    
                    finally:
                        field.blockSignals(False)
    
                # ----------------------------
                # QComboBox
                # ----------------------------
                elif isinstance(field, QComboBox):
    
                    field.blockSignals(True)
    
                    try:
    
                        if value is None:
    
                            field.setCurrentIndex(-1)
    
                        else:
    
                            idx = field.findText(str(value))
    
                            if idx >= 0:
                                field.setCurrentIndex(idx)
    
                    finally:
    
                        field.blockSignals(False)
    
            # Remove extra rows
            while len(self.rows) > len(values):
    
                row_widget, _, _ = self.rows[-1]
    
                self.remove_row(row_widget)
    
        finally:
    
            self._internal_update = False
    
        if not self._internal_update:
            self.value_changed.emit()

# ----------------------------
# Dependency Widget
# ----------------------------

class DependencyWidget(QWidget, BaseParamWidget):
    value_changed = Signal()

    def __init__(self,
             parent_widget,
             child_widget,
             metadata_dict,
             child_source,
             use_dependency=True,
             allow_empty=False):

        super().__init__()

        self.parent_widget = parent_widget
        self.child_widget = child_widget
        self.metadata_dict = metadata_dict
        self.child_source = child_source if isinstance(child_source, list) else [child_source]
        self.allow_empty = allow_empty
        self.use_dependency = use_dependency

        # initial fill
        self.update_child()

        # dynamic update
        self.parent_widget.currentTextChanged.connect(self.update_child)
        self.parent_widget.currentTextChanged.connect(
            lambda: self.value_changed.emit()
        )
        
        self._apply_default_size(self.parent_widget)
        
        # Child signal
        if isinstance(self.child_widget, QComboBox):
            self.child_widget.currentTextChanged.connect(
                lambda: self.value_changed.emit()
            )
        
        elif isinstance(self.child_widget, QLineEdit):
            self.child_widget.textChanged.connect(
                lambda: self.value_changed.emit()
            )
        
        elif hasattr(self.child_widget, "value_changed"):
            self.child_widget.value_changed.connect(
                lambda: self.value_changed.emit()
            )

    def update_child(self):
        if not self.use_dependency:
            return

        parent_value = self.parent_widget.currentText()
    
        all_values = []
    
        for src in self.child_source:
    
            if src == "__dummy_type__":
                values = self.metadata_dict.get("__dummy_type__", [])
            else:
                values = self.metadata_dict.get(parent_value, [])
    
            all_values.extend(values)
    
        seen = set()
        unique_values = []
        for v in all_values:
            if v not in seen:
                unique_values.append(v)
                seen.add(v)
    
        if isinstance(self.child_widget, QComboBox):
            self.child_widget.clear()
            self.child_widget.addItems([str(v) for v in unique_values])
            self.child_widget.setCurrentIndex(-1)
    
        elif isinstance(self.child_widget, QLineEdit):
            self.child_widget.setText(str(unique_values[0]) if unique_values else "")
            self.child_widget.setMinimumWidth(140)
            self.child_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
    
        elif hasattr(self.child_widget, "set_value"):
            self.child_widget.set_value(unique_values)
            
        self._apply_default_size(self.child_widget)
    
        QTimer.singleShot(0, self.value_changed.emit)
            
    def get_value(self):

        def read(widget):
            if hasattr(widget, "get_value"):
                return widget.get_value()
            if isinstance(widget, QComboBox):
                return widget.currentText()
            if isinstance(widget, QLineEdit):
                return widget.text()
            return None
    
        return {
            "parent": read(self.parent_widget),
            "child": read(self.child_widget),
            "child_source": self.child_source
        }
    
    def set_value(self, value):

        if not isinstance(value, dict):
            return
    
        parent = value.get("parent")
        child = value.get("child")
    
        # set parent first
        self.parent_widget.set_value(parent)
    
        # refresh child entries
        self.update_child()
    
        # then child
        self.child_widget.set_value(child)
    
    
    def clear_value(self):
    
        self.parent_widget.clear_value()
        self.update_child()
        self.child_widget.clear_value()
        
    def _apply_default_size(self, widget):
        widget.setMinimumWidth(140)
        widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)

# ----------------------------
# Dependency RowWidget (subwidget)
# ----------------------------

class DependencyRowWidget(QWidget, BaseParamWidget):
    value_changed = Signal() 

    def __init__(self,
             metadata_dict,
             data_list,
             child_metas,
             parent_type,
             create_widget_func,
             parse_gui_type_func,
             allow_empty=False):

        super().__init__()
        self._internal_update = False

        self.metadata_dict = metadata_dict
        self.data_list = data_list
        self.child_metas = child_metas
        self.parent_type = parent_type
        self.create_widget_func = create_widget_func
        self.parse_gui_type_func = parse_gui_type_func
        self.allow_empty = allow_empty

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0,0,0,0)

        # parent
        self.parent_combo = ComboBoxAdapter()

        if self.allow_empty:
            self.parent_combo.addItem("")
                
        if self.parent_type[0] == "meta":
            keys = [
                k for k in metadata_dict.keys()
                if k != "__dummy_type__"
            ]
        elif self.parent_type[0] == "data":
            keys = [
                k for k in data_list
            ]
    
        self.parent_combo.addItems(keys)
        self.parent_combo.setCurrentIndex(-1)
        
        self.parent_combo.currentTextChanged.connect(
            lambda *_: self.value_changed.emit()
        )
        self.parent_combo.currentTextChanged.connect(
            lambda *_: self.update_child()
        )
        
        self._apply_default_size(self.parent_combo)

        layout.addWidget(self.parent_combo)

        self.child_widgets = []
        
        for child_meta in self.child_metas:
        
            child_type = self.parse_gui_type_func(
                child_meta["GUI"],
                child_meta.get("Options")
            )
        
            child_widget = self.create_widget_func(
                child_type,
                child_meta["GUI"]
            )
            
            self._apply_default_size(child_widget)
        
            layout.addWidget(child_widget)
        
            self.child_widgets.append(
                (
                    child_widget,
                    child_meta
                )
            )
        
            # signals
            if isinstance(child_widget, QComboBox):
                child_widget.currentTextChanged.connect(
                    lambda *_: self.value_changed.emit()
                )
        
            elif isinstance(child_widget, QLineEdit):
                child_widget.textChanged.connect(
                    lambda *_: self.value_changed.emit()
                )
        
            elif hasattr(child_widget, "value_changed"):
                child_widget.value_changed.connect(
                    lambda: self.value_changed.emit()
                )

        # initial fill
        self.update_child()

    def update_child(self):
    
        if self._internal_update:
            return
    
        self._internal_update = True
    
        try:
            parent_column = self.parent_combo.get_value()
    
            for child_widget, child_meta in self.child_widgets:
    
                gui = child_meta.get("GUI", "")
                use_dependency = "X" in gui
    
                child_source = child_meta.get("M_column")
    
                sources = (
                    child_source
                    if isinstance(child_source, list)
                    else [child_source]
                )
    

                all_values = []
    
                for src in sources:
    
                    if use_dependency:
                        # dependency mode -> always use parent mapping
                        if src == "__dummy_type__":
                            values = self.metadata_dict.get("__dummy_type__", [])
                        else:
                            values = self.metadata_dict.get(parent_column, [])
    
                    else:
                        # non-dependency mode -> direct source mapping
                        values = self.metadata_dict.get(src, [])
    
                    all_values.extend(values)
    
                # dedupe
                seen = set()
                unique_values = []
    
                for v in all_values:
                    if v not in seen:
                        unique_values.append(v)
                        seen.add(v)
    
                # ----------------------------
                # COMBOBOX
                # ----------------------------
                if isinstance(child_widget, QComboBox):
    
                    if use_dependency:
                        current = child_widget.currentText()
    
                        child_widget.blockSignals(True)
                        child_widget.clear()
                        child_widget.addItems([str(v) for v in unique_values])
    
                        idx = child_widget.findText(current)
                        child_widget.setCurrentIndex(idx if idx >= 0 else -1)
    
                        child_widget.blockSignals(False)
    
                    else:
                        pass
    
                # ----------------------------
                # LINE EDIT
                # ----------------------------
                elif isinstance(child_widget, QLineEdit):
    
                    if use_dependency:
                        new_text = str(unique_values[0]) if unique_values else ""
    
                        if child_widget.text() != new_text:
                            child_widget.blockSignals(True)
                            child_widget.setText(new_text)
                            child_widget.blockSignals(False)
    
                    else:
                        pass
    
                # ----------------------------
                # MULTI VALUE WIDGET
                # ----------------------------
                elif isinstance(child_widget, MultiValueWidget):
    
                    if not use_dependency:
                        continue
    
                    for row_widget, _, _ in child_widget.rows:
    
                        inner = row_widget.layout().itemAt(0).widget()
    
                        # ComboBox rows
                        if isinstance(inner, QComboBox):
    
                            current = inner.currentText()
    
                            inner.blockSignals(True)
                            inner.clear()
                            inner.addItems([str(v) for v in unique_values])
    
                            idx = inner.findText(current)
                            inner.setCurrentIndex(idx if idx >= 0 else -1)
    
                            inner.blockSignals(False)
    
                        # LineEdit rows
                        elif isinstance(inner, QLineEdit):
    
                            new_text = str(unique_values[0]) if unique_values else ""
    
                            if inner.text() != new_text:
                                inner.blockSignals(True)
                                inner.setText(new_text)
                                inner.blockSignals(False)
    
        finally:
            self._internal_update = False
    
        QTimer.singleShot(0, self.value_changed.emit)
            
    def get_value(self):

        result = {
            "parent": self.parent_combo.get_value(),
            "children": []
        }
    
        for child_widget, child_meta in self.child_widgets:
    
            if hasattr(child_widget, "get_value"):
                child_val = child_widget.get_value()
    
            elif isinstance(child_widget, QComboBox):
                child_val = child_widget.currentText()
    
            elif isinstance(child_widget, QLineEdit):
                child_val = child_widget.text()
    
            else:
                child_val = None
    
            result["children"].append({
                "child": child_val,
                "child_source": child_meta.get("M_column")
            })
    
        return result
        
    def set_value(self, value):
    
        if self._internal_update:
            return
    
        if not isinstance(value, dict):
            return
    
        self._internal_update = True
    
        try:
    
            self.parent_combo.blockSignals(True)
            self.parent_combo.set_value(value.get("parent"))
            self.parent_combo.blockSignals(False)
    
            self.update_child()
    
            children = value.get("children", [])
    
            for i, child_data in enumerate(children):
    
                if i >= len(self.child_widgets):
                    break
    
                child_widget, _ = self.child_widgets[i]
    
                child_val = child_data.get("child")
    
                if hasattr(child_widget, "set_value"):
    
                    child_widget.blockSignals(True)
                    child_widget.set_value(child_val)
                    child_widget.blockSignals(False)
    
                elif isinstance(child_widget, QComboBox):
    
                    child_widget.blockSignals(True)
                    child_widget.setCurrentText(str(child_val))
                    child_widget.blockSignals(False)
    
                elif isinstance(child_widget, QLineEdit):
    
                    child_widget.blockSignals(True)
                    child_widget.setText(
                        "" if child_val is None else str(child_val)
                    )
                    child_widget.blockSignals(False)
    
        finally:
    
            self._internal_update = False
    
        self.value_changed.emit()
    
    def clear_value(self):
    
        self.parent_combo.clear_value()
    
        self.update_child()
    
        for child_widget, _ in self.child_widgets:
    
            if hasattr(child_widget, "clear_value"):
    
                child_widget.clear_value()
                
    def _apply_default_size(self, widget):
        widget.setMinimumWidth(140)
        widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            
# ----------------------------
# MultiDependency Widget
# ----------------------------
        
class MultiDependencyWidget(QWidget, BaseParamWidget):

    value_changed = Signal()

    def __init__(self,
             metadata_dict,
             data_list,
             parent_type,
             child_metas,
             create_widget_func,
             parse_gui_type_func,
             allow_empty=False):

        super().__init__()
        self._internal_update = False

        self.metadata_dict = metadata_dict
        self.data_list = data_list
        self.parent_type = parent_type
        self.child_metas = child_metas
        self.create_widget_func = create_widget_func
        self.parse_gui_type_func = parse_gui_type_func
        self.allow_empty = allow_empty

        self.layout = QVBoxLayout(self)

        self.layout.setContentsMargins(0,0,0,0)
        self.layout.setSpacing(2)

        self.rows = []

        self.add_row()

    def add_row(self):

        row_widget = QWidget()

        row_layout = QHBoxLayout(row_widget)

        row_layout.setContentsMargins(0,0,0,0)

        dep_widget = DependencyRowWidget(
            metadata_dict=self.metadata_dict,
            data_list=self.data_list,
            child_metas=self.child_metas,
            parent_type=self.parent_type,
            create_widget_func=self.create_widget_func,
            parse_gui_type_func=self.parse_gui_type_func,
            allow_empty=self.allow_empty
        )
        
        dep_widget.value_changed.connect(
            lambda: (
                None if self._internal_update
                else self.value_changed.emit()
            )
        )

        add_btn = QPushButton("+")
        remove_btn = QPushButton("-")

        add_btn.setFixedWidth(25)
        remove_btn.setFixedWidth(25)

        add_btn.clicked.connect(self.add_row)

        remove_btn.clicked.connect(
            lambda: self.remove_row(row_widget)
        )

        row_layout.addWidget(dep_widget, 1)
        row_layout.addWidget(add_btn)
        row_layout.addWidget(remove_btn)

        self.layout.addWidget(row_widget)

        self.rows.append(
            (row_widget, add_btn, remove_btn)
        )

        self.update_buttons()

    def remove_row(self, row_widget):

        if len(self.rows) <= 1:
            return
    
        for i, (row, _, _) in enumerate(self.rows):
    
            if row == row_widget:
    
                self.rows.pop(i)
    
                row.setParent(None)
                row.deleteLater()
    
                break
    
        self.update_buttons()
    
        if not self._internal_update:
            self.value_changed.emit()

    def update_buttons(self):

        n = len(self.rows)

        for i, (row, add_btn, remove_btn) in enumerate(self.rows):

            add_btn.setEnabled(i == n - 1)

            remove_btn.setVisible(n > 1)
            
    def is_valid(self):

        for row_widget, _, _ in self.rows:
            dep_widget = row_widget.layout().itemAt(0).widget()
            parent_val = dep_widget.parent_combo.get_value()
            if not parent_val:
                return False
            for child_widget, _ in dep_widget.child_widgets:
                if hasattr(child_widget, "get_value"):
                    child_val = child_widget.get_value()
                elif isinstance(child_widget, QComboBox):
                    child_val = child_widget.currentText()
                elif isinstance(child_widget, QLineEdit):
                    child_val = child_widget.text()
                else:
                    child_val = None
                if not child_val:
                    return False
    
        return True
            
    def get_value(self):

        values = []
    
        for row_widget, _, _ in self.rows:
            dep_widget = row_widget.layout().itemAt(0).widget()
            children = []
            for child_widget, child_meta in dep_widget.child_widgets:
                if hasattr(child_widget, "get_value"):
                    child_val = child_widget.get_value()
                elif isinstance(child_widget, QComboBox):
                    child_val = child_widget.currentText()
                elif isinstance(child_widget, QLineEdit):
                    child_val = child_widget.text()
                else:
                    child_val = None
    
                children.append({
                    "child": child_val,
                    "child_source": child_meta.get("M_column")
                })
    
            values.append({
                "parent": dep_widget.parent_combo.get_value(),
                "children": children
            })
    
        return values
    
    def set_value(self, values):
        
        if getattr(self, "_syncing", False):
            return

        if self._internal_update:
            return
    
        self._internal_update = True
        try:
    
            if values is None:
                values = []
    
            if not isinstance(values, list):
                values = [values]
    
            while len(self.rows) < len(values):
                self.add_row()
    
            for i, value in enumerate(values):
    
                row_widget, _, _ = self.rows[i]
                dep_widget = row_widget.layout().itemAt(0).widget()
    
                if isinstance(value, dict) and "child_sources" in value:
                    dep_widget.child_source = value["child_sources"]
    
                dep_widget.update_child()
    
                dep_widget.set_value(value)
    
            while len(self.rows) > max(1, len(values)):
                row_widget, _, _ = self.rows[-1]
                self.remove_row(row_widget)
    
        finally:
            self._internal_update = False
    
    def clear_value(self):
        self.set_value([])

# ----------------------------
# Config: General
# ----------------------------

class GeneralParameterConfig(QWidget):
    def __init__(self, step, flags_dict, params_dict, context):
        super().__init__()

        layout = QVBoxLayout(self)

        widget = DynamicParameterConfig(
            step,
            flags_dict,
            params_dict,
            meta_columns=context.get("metadata", {}),
            data_columns=context.get("data_header", [])
        )

        layout.addWidget(widget)

# ----------------------------
# Config: Filtering_long
# ----------------------------

class FilteringLongParameterConfig(QWidget):
    def __init__(self, step, flags_dict, params_dict, context):
        super().__init__()

        layout = QVBoxLayout(self)
        
        flags_dict["FILTERING_LONG"]["input_origin"]["PriorityFirst"] = True
        flags_dict["FILTERING_LONG"]["filter_proteotypic"]["VisibleIf"] = {"input_origin":"diann"}
        flags_dict["FILTERING_LONG"]["filter_by_qvalue"]["VisibleIf"] = {"input_origin":"diann"}
        flags_dict["FILTERING_LONG"]["filter_by_global_qvalue"]["VisibleIf"] = {"input_origin":"diann"}
        flags_dict["FILTERING_LONG"]["filter_by_gg_qvalue"]["VisibleIf"] = {"input_origin":"diann"}
        flags_dict["FILTERING_LONG"]["filter_by_pg_qvalue"]["VisibleIf"] = {"input_origin":"diann"}
        flags_dict["FILTERING_LONG"]["filter_by_global_pg_qvalue"]["VisibleIf"] = {"input_origin":"diann"}
        flags_dict["FILTERING_LONG"]["filter_by_p_qvalue"]["VisibleIf"] = {"input_origin":"diann"}
        flags_dict["FILTERING_LONG"]["filter_by_lib_qvalue"]["VisibleIf"] = {"input_origin":"diann"}
        flags_dict["FILTERING_LONG"]["filter_by_lib_pg_qvalue"]["VisibleIf"] = {"input_origin":"diann"}
        flags_dict["FILTERING_LONG"]["filter_by_PEP"]["VisibleIf"] = {"input_origin":"diann"}
        flags_dict["FILTERING_LONG"]["filter_by_MScorr"]["VisibleIf"] = {"input_origin":"diann"}
        flags_dict["FILTERING_LONG"]["filter_by_zero_values"]["VisibleIf"] = {"input_origin":"diann"}
        flags_dict["FILTERING_LONG"]["filter_by_evidence"]["VisibleIf"] = {"input_origin":"diann"}
        flags_dict["FILTERING_LONG"]["filter_by_MScorr_evidence_combined"]["VisibleIf"] = {"input_origin":"diann"}
        flags_dict["FILTERING_LONG"]["filter_by_precursorcharge"]["VisibleIf"] = {"input_origin":"diann"}

        widget = DynamicParameterConfig(
            step,
            flags_dict,
            params_dict,
            meta_columns=context.get("metadata", {}),
            data_columns=context.get("data_header", [])
        )

        layout.addWidget(widget)

# ----------------------------
# Config: Filtering wide
# ----------------------------

class FilteringWideParameterConfig(QWidget):
    def __init__(self, step, flags_dict, params_dict, context):
        super().__init__()

        layout = QVBoxLayout(self)

        widget = DynamicParameterConfig(
            step,
            flags_dict,
            params_dict,
            meta_columns=context.get("metadata", {}),
            data_columns=context.get("data_header", [])
        )

        layout.addWidget(widget)
        
# ----------------------------
# Config: Normalisation
# ----------------------------

class NormalisationParameterConfig(QWidget):
    def __init__(self, step, flags_dict, params_dict, context):
        super().__init__()

        layout = QVBoxLayout(self)

        flags_dict["NORMALISATION"]["on_peptide_list"]["VisibleIf"] = {"normalization_method":"median_scaling"}
        params_dict["normalization_method"]["median_scaling"]["nan_intersect_thresh"]["VisibleIf"] = {"normalization_on_intersecting_precs":True}
        params_dict["normalization_method"]["median_scaling"]["on_mad"]["VisibleIf"] = {"exclude_outlier_nfs":True}
        params_dict["normalization_method"]["median_scaling"]["rel_threshold"]["VisibleIf"] = {"on_mad":True}
        params_dict["normalization_method"]["median_scaling"]["min_nf"]["VisibleIf"] = {"on_mad":False,
                                                                                        "exclude_outlier_nfs":True}
        params_dict["normalization_method"]["median_scaling"]["max_nf"]["VisibleIf"] = {"on_mad":False,
                                                                                        "exclude_outlier_nfs":True}

        widget = DynamicParameterConfig(
            step,
            flags_dict,
            params_dict,
            meta_columns=context.get("metadata", {}),
            data_columns=context.get("data_header", [])
        )

        layout.addWidget(widget)
        
# ----------------------------
# Config: Imputation
# ----------------------------
        
class ImputationParameterConfig(QWidget):
    def __init__(self, step, flags_dict, params_dict, context):
        super().__init__()
        
        layout = QVBoxLayout(self)
        # Move to position one
        flags_dict["IMPUTATION"]["simple_imputation"]["PriorityFirst"] = True
        # Exclude flags
        flags_dict["IMPUTATION"]["imputation_method_HA"]["VisibleIf"] = {"simple_imputation":False}
        flags_dict["IMPUTATION"]["imputation_method_LA"]["VisibleIf"] = {"simple_imputation":False}
        flags_dict["IMPUTATION"]["imputation_method_ProtAbsence"]["VisibleIf"] = {"simple_imputation":False}
        flags_dict["IMPUTATION"]["imputation_method_ProtAbsenceHA"]["VisibleIf"] = {"simple_imputation":False}
        flags_dict["IMPUTATION"]["imputation_method_ProtAbsenceLA"]["VisibleIf"] = {"simple_imputation":False}
        flags_dict["IMPUTATION"]["eval_missingness_on_protein_presence"]["VisibleIf"] = {"simple_imputation":False}
        flags_dict["IMPUTATION"]["eval_missingness_by_detection_threshold"]["VisibleIf"] = {"simple_imputation":False}
        flags_dict["IMPUTATION"]["hierarchical_imputation"]["VisibleIf"] = {"simple_imputation":False}
        flags_dict["IMPUTATION"]["remember_missingness_type_in_mask"]["VisibleIf"] = {"simple_imputation":False}
        flags_dict["IMPUTATION"]["save_intermediate_output_for_QA"]["VisibleIf"] = {"simple_imputation":False}
        flags_dict["IMPUTATION"]["impute_over_whole_dataset"]["VisibleIf"] = {"simple_imputation":False}
        # Exclude values
        flags_dict["IMPUTATION"]["imputation_method"]["FilteredOptionsIf"] = {"simple_imputation":{False:["missforest","bpca_R","ppca","ppm","lcg"]}}
        flags_dict["IMPUTATION"]["imputation_method_HA"]["FilteredOptionsIf"] = {"simple_imputation":{False:["missforest","bpca_R","ppca","ppm","lcg"]}}
        flags_dict["IMPUTATION"]["imputation_method_LA"]["FilteredOptionsIf"] = {"simple_imputation":{False:["missforest","bpca_R","ppca","ppm","lcg"]}}
        flags_dict["IMPUTATION"]["imputation_method_ProtAbsence"]["FilteredOptionsIf"] = {"simple_imputation":{False:["missforest","bpca_R","ppca","ppm","lcg"]}}
        flags_dict["IMPUTATION"]["imputation_method_ProtAbsenceHA"]["FilteredOptionsIf"] = {"simple_imputation":{False:["missforest","bpca_R","ppca","ppm","lcg"]}}
        flags_dict["IMPUTATION"]["imputation_method_ProtAbsenceLA"]["FilteredOptionsIf"] = {"simple_imputation":{False:["missforest","bpca_R","ppca","ppm","lcg"]}}
                
        widget = DynamicParameterConfig(
            step,
            flags_dict,
            params_dict,
            meta_columns=context.get("metadata", {}),
            data_columns=context.get("data_header", [])
        )

        layout.addWidget(widget)
        
# ----------------------------
# Config: Correction
# ----------------------------
        
class CorrectionParameterConfig(QWidget):
    def __init__(self, step, flags_dict, params_dict, context):
        super().__init__()

        layout = QVBoxLayout(self)
        
        flags_dict["CORRECTION"]["check_feature_presence_batch"]["VisibleIf"] = {"correct_method":["combat","harmony"]}
        
        widget = DynamicParameterConfig(
            step,
            flags_dict,
            params_dict,
            meta_columns=context.get("metadata", {}),
            data_columns=context.get("data_header", [])
        )

        layout.addWidget(widget)
        
# ----------------------------
# Config: Summarisation
# ----------------------------
        
class SummarisationParameterConfig(QWidget):
    def __init__(self, step, flags_dict, params_dict, context):
        super().__init__()

        layout = QVBoxLayout(self)

        widget = DynamicParameterConfig(
            step,
            flags_dict,
            params_dict,
            meta_columns=context.get("metadata", {}),
            data_columns=context.get("data_header", [])
        )

        layout.addWidget(widget)

# ----------------------------
# Config Factory
# ----------------------------

class StepConfigFactory:
    registry = {}

    @classmethod
    def register(cls, name, widget):
        cls.registry[name] = widget

    @classmethod
    def create(cls, step, flags_dict, params_dict, context):
        widget_class = cls.registry.get(step.name)

        if widget_class:
            return widget_class(step, flags_dict, params_dict, context)
        return QLabel("No configuration available")

StepConfigFactory.register("General", GeneralParameterConfig)
StepConfigFactory.register("Filtering_long", FilteringLongParameterConfig)
StepConfigFactory.register("Filtering_wide", FilteringWideParameterConfig)
StepConfigFactory.register("Normalisation", NormalisationParameterConfig)
StepConfigFactory.register("Imputation", ImputationParameterConfig)
StepConfigFactory.register("Correction", CorrectionParameterConfig)
StepConfigFactory.register("Summarisation", SummarisationParameterConfig)

# ----------------------------
# Line Edit Adapter
# ----------------------------

class LineEditAdapter(QLineEdit, BaseParamWidget):

    def clear_value(self):
        self.clear()

    def set_value(self, value):
        self.setText("" if value is None else str(value))

    def get_value(self):
        text = self.text().strip()
        return text if text else None


# ----------------------------
# Combo Box Adapter
# ----------------------------

class ComboBoxAdapter(QComboBox, BaseParamWidget):

    def clear_value(self):
        self.setCurrentIndex(-1)

    def set_value(self, value):

        if value is None:
            self.setCurrentIndex(-1)
            return

        idx = self.findText(str(value))
        self.setCurrentIndex(idx if idx >= 0 else -1)

    def get_value(self):

        if self.currentIndex() < 0:
            return None

        text = self.currentText()

        if text in ("", "None"):
            return None

        return text


# ----------------------------
# Checkbox Adapter
# ----------------------------

class CheckBoxAdapter(QCheckBox, BaseParamWidget):

    def clear_value(self):
        self.setChecked(False)

    def set_value(self, value):
        # Normalize bool-like values
        if isinstance(value, str):
            value = value.strip().lower()
            if value in ("1",1):
                checked = True
            elif value in ("0",0):
                checked = False
            else:
                checked = bool(value)
        else:
            checked = bool(value)

        self.setChecked(checked)

    def get_value(self):
        return str(int(self.isChecked()))

# ----------------------------
# File Browser Adapter
# ----------------------------
    
class FileBrowserAdapter(QWidget, BaseParamWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.line = QLineEdit(self)
        self.button = QPushButton("Browse")

        layout = QHBoxLayout(self)
        layout.addWidget(self.line)
        layout.addWidget(self.button)

        self.button.clicked.connect(self.open_dialog)

    def open_dialog(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select File")
        if path:
            self.line.setText(path)

    def clear_value(self):
        self.line.clear()

    def set_value(self, value):
        self.line.blockSignals(True)
        self.line.setText("" if value is None else str(value))
        self.line.blockSignals(False)

    def get_value(self):
        text = self.line.text().strip()
        return text if text else None
    
# ----------------------------
# Range Widget Adapter
# ----------------------------

class RangeWidgetAdapter(QWidget, BaseParamWidget):

    def __init__(self, integer=False):
        super().__init__()

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(2)

        self.min_field = LineEditAdapter()
        self.max_field = LineEditAdapter()
        self.integer = integer

        if self.integer:
            self.min_field.setValidator(QIntValidator())
            self.max_field.setValidator(QIntValidator())
        else:
            self.min_field.setValidator(QDoubleValidator())
            self.max_field.setValidator(QDoubleValidator())

        layout.addWidget(self.min_field)
        layout.addWidget(self.max_field)

    def clear_value(self):
        self.min_field.clear()
        self.max_field.clear()

    def set_value(self, value):

        if value is None:
            self.clear_value()
            return

        min_val, max_val = value

        self.min_field.set_value(min_val)

        if max_val == np.inf:
            self.max_field.clear()
        else:
            self.max_field.set_value(max_val)

    def get_value(self):

        min_text = self.min_field.text().strip()
        max_text = self.max_field.text().strip()
        
        def _to_float(text):
            try:
                return float(text)
            except ValueError:
                return None
        
        
        if self.integer:
            min_val = int(min_text) if min_text not in ("", "-", "+") else None
            max_val = int(max_text) if max_text not in ("", "-", "+") else None
        else:
            min_val = _to_float(min_text)
            max_val = _to_float(max_text)

        if min_val is None and max_val is None:
            return None

        if max_val is None:
            max_val = np.inf

        return (min_val, max_val)
    
    def is_valid(self):
        min_text = self.min_field.text().strip()
        max_text = self.max_field.text().strip()
    
        def ok(x):
            try:
                float(x)
                return True
            except:
                return False
    
        if min_text in ("", "-", "+") or max_text in ("", "-", "+"):
            return False
    
        if min_text and not ok(min_text):
            return False
        if max_text and not ok(max_text):
            return False
    
        if ok(min_text) and ok(max_text):
            if float(max_text) < float(min_text):
                return False
    
        return True
    
class TripleIntWidget(QWidget, BaseParamWidget):

    def __init__(self):
        super().__init__()

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0,0,0,0)

        self.f1 = LineEditAdapter()
        self.f2 = LineEditAdapter()
        self.f3 = LineEditAdapter()

        validator = QIntValidator()

        self.f1.setValidator(validator)
        self.f2.setValidator(validator)
        self.f3.setValidator(validator)

        layout.addWidget(self.f1)
        layout.addWidget(self.f2)
        layout.addWidget(self.f3)

    def get_value(self):
        vals = []

        for f in (self.f1, self.f2, self.f3):
            t = f.text().strip()
            vals.append(int(t) if t else None)

        return tuple(vals)

    def set_value(self, value):

        if not value:
            return

        for field, val in zip((self.f1, self.f2, self.f3), value):

            field.setText("" if val is None else str(val))

    def clear_value(self):

        self.f1.clear()
        self.f2.clear()
        self.f3.clear()
        
    def is_valid(self):
        for f in (self.f1, self.f2, self.f3):
            t = f.text().strip()
            if not t:
                return False
            try:
                int(t)
            except:
                return False
        return True

# ----------------------------
# Config Panel
# ----------------------------
class ConfigPanel(QWidget):
    def __init__(self, flags_dict, params_dict, context):
        super().__init__()
        self.context = context
        self.flags_dict = flags_dict
        self.params_dict = params_dict

        outer_layout = QVBoxLayout(self)
        outer_layout.setContentsMargins(0, 0, 0, 0)
        outer_layout.setSpacing(0)

        self.scroll = QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.scroll.setFrameShape(QFrame.NoFrame)

        self.container = QWidget()
        self.layout = QVBoxLayout(self.container)
        self.layout.setContentsMargins(4, 4, 4, 4)
        self.layout.setSpacing(4)
        self.layout.setAlignment(Qt.AlignTop)

        self.scroll.setWidget(self.container)
        outer_layout.addWidget(self.scroll)

        self.layout.addWidget(QLabel("Select a step to configure"))

    def load_config_for_step(self, current):
        # Clear layout
        for i in reversed(range(self.layout.count())):
            widget = self.layout.itemAt(i).widget()
            if widget:
                widget.deleteLater()
    
        if current:
            step = current.data(Qt.UserRole)
    
            label = QLabel(f"Config: {step.name}")
            label.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)
    
            widget = StepConfigFactory.create(step, self.flags_dict, self.params_dict, self.context)
            widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Maximum)
    
            self.layout.addWidget(label)
            self.layout.addWidget(widget)
    
        else:
            self.layout.addWidget(QLabel("Select a step to configure"))


# ----------------------------
# Global Settings
# ----------------------------

class GlobalSettingsPanel(QWidget):
    settings_valid = Signal(bool)

    def __init__(self, context):
        super().__init__()

        self.context = context
        self.layout = QVBoxLayout(self)
        self.layout.setSpacing(2)
        self.layout.setContentsMargins(2, 2, 2, 2)

        # ----------------------------
        # Create fields
        # ----------------------------

        #### For testing and debugging
        #self.ms_file = QLineEdit("./data/data_00.tsv")
        #self.meta_file = QLineEdit("./data/Metadata.tsv")
        #self.param_dir = QLineEdit("../")
        #self.work_dir = QLineEdit("./")
        
        self.ms_file = QLineEdit("")
        self.meta_file = QLineEdit("")
        self.param_dir = QLineEdit("")
        self.work_dir = QLineEdit("")

        self.ms_file.setPlaceholderText("Select proteomic data file...")
        self.meta_file.setPlaceholderText("Select metadata file...")
        self.param_dir.setPlaceholderText("Select directory with parameter dictionaries...")
        self.work_dir.setPlaceholderText("Select working directory...")

        self.config_file = QLineEdit()
        self.config_file.setPlaceholderText("Select old configuration file...")
        self.config_file.setEnabled(False)

        self.dropdown = QComboBox()
        self.dropdown.setEnabled(False)

        self.rerun_checkbox = QCheckBox("Rerun?")
        self.docker_checkbox = QCheckBox("Docker")
        self.docker_checkbox.setChecked(True)

        self.config_file_button = QPushButton("Browse")
        self.config_file_button.setEnabled(False)
        self.config_file_button.clicked.connect(self.select_config_file)

        # ----------------------------
        # Helper
        # ----------------------------
        def make_file_row(line_edit, button, callback=None):
            if callback:
                button.clicked.connect(lambda _, le=line_edit: callback(le))
            row = QHBoxLayout()
            row.addWidget(line_edit)
            row.addWidget(button)
            row.setSpacing(2)
            row.setContentsMargins(0, 0, 0, 0)
            return row

        input_box = QGroupBox("Input / Paths")
        grid = QGridLayout()

        grid.addWidget(QLabel("MS-Data"), 0, 0)
        grid.addLayout(make_file_row(self.ms_file, QPushButton("Browse"), self.select_file), 1, 0)

        grid.addWidget(QLabel("Metadata"), 2, 0)
        grid.addLayout(make_file_row(self.meta_file, QPushButton("Browse"), self.select_file), 3, 0)

        grid.addWidget(QLabel("Parameter dictionaries"), 0, 1)
        grid.addLayout(make_file_row(self.param_dir, QPushButton("Browse"), self.select_dir), 1, 1)

        grid.addWidget(QLabel("Working directory"), 2, 1)
        grid.addLayout(make_file_row(self.work_dir, QPushButton("Browse"), self.select_dir), 3, 1)

        input_box.setLayout(grid)

        rerun_box = QGroupBox("Rerun options + Docker")
        rerun_layout = QVBoxLayout()

        checkbox_layout = QHBoxLayout()
        checkbox_layout.addWidget(self.rerun_checkbox)
        checkbox_layout.addWidget(self.docker_checkbox)
        checkbox_layout.insertStretch(1)
        
        rerun_layout.addLayout(checkbox_layout)
        rerun_layout.addWidget(QLabel("Old configuration file"))
        rerun_layout.addLayout(make_file_row(self.config_file, self.config_file_button, None))
        rerun_layout.addWidget(self.dropdown)

        grid.setSpacing(2)
        grid.setContentsMargins(2, 2, 2, 2)
        
        rerun_layout.setSpacing(2)
        rerun_layout.setContentsMargins(2, 2, 2, 2)

        rerun_box.setLayout(rerun_layout)

        self.rerun_checkbox.stateChanged.connect(self.toggle_rerun_fields)

        self.load_button = QPushButton("Load Data")
        self.load_button.clicked.connect(self.validate)

        self.status_label = QLabel("Not loaded")

        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.setVisible(False)

        top_layout = QHBoxLayout()
        top_layout.addWidget(input_box, 4)
        top_layout.addWidget(rerun_box, 1)
        top_layout.setSpacing(10)
        
        self.layout.addLayout(top_layout)
        self.layout.addWidget(self.load_button)
        self.layout.addWidget(self.status_label)
        self.layout.addWidget(self.progress)
        self.layout.addWidget(self.load_button)
        self.layout.addWidget(self.status_label)
        self.layout.addWidget(self.progress)

    def select_file(self, target):
        path, _ = QFileDialog.getOpenFileName(self, dir=os.getcwd())
        if path:
            target.setText(path)

    def select_dir(self, target):
        path = QFileDialog.getExistingDirectory(self, dir=os.getcwd())
        if path:
            target.setText(path)

    def validate(self):
        valid = all([
            self.ms_file.text(),
            #self.meta_file.text(), ### Can be empty
            self.param_dir.text(),
            self.work_dir.text()
        ])

        if not valid:
            self.status_label.setText("Missing inputs ❌")
            self.settings_valid.emit(False)
            return

        self.status_label.setText("Loading...")
        self.progress.setVisible(True)

        QTimer.singleShot(100, self._do_loading)

    def _do_loading(self):
        try:
            self.load_data()
            if not self.flag_meta:
                self.status_label.setText("Failed loading metadata ❌")
            elif not self.flag_ms:
                self.status_label.setText("Failed loading MS data ❌")
            elif not self.flag_params:
                self.status_label.setText("Failed loading parameters ❌")
            else:
                self.status_label.setText("Loaded ✅")
                self.settings_valid.emit(True)
                main = self.window()

                if hasattr(main, "workflow_panel") and hasattr(main, "config_panel"):
                
                    current = main.workflow_panel.workflow_list.currentItem()
                
                    if current is not None:
                        main.config_panel.load_config_for_step(current)
        except Exception as e:
            #import traceback
            #traceback.print_exc()
            print(e)
            self.status_label.setText("Failed ❌")
            self.settings_valid.emit(False)

        self.progress.setVisible(False)

    def load_data(self):
        try:
            self.metadata = gui_funcs.load_metadata(self.meta_file.text())
            self.flag_meta = self.metadata is not None
            
            if self.flag_meta:
                meta_dict = {}
            
                for col in self.metadata.columns:
                    # unique values, drop NaN
                    values = self.metadata[col].dropna().unique().tolist()
                    meta_dict[col] = values
            
                meta_dict["__dummy_type__"] = ["ordinal", "nominal", "continuous"]
            
                # store full dict
                self.context["metadata"] = meta_dict
            else:
                self.context["metadata"] = {"File.Name":[]}
                if not self.meta_file.text():
                    self.flag_meta = True
    
        except:
            self.flag_meta = False
            self.context["metadata"] = []
    
        try:
            self.ms_header = gui_funcs.load_data_header(self.ms_file.text())
            self.flag_ms = True

            self.context["data_header"] = self.ms_header or []

        except:
            self.flag_ms = False
            self.context["data_header"] = []
    
        try:
            self.flags_dict, self.params_dict = gui_funcs.load_parameter_dictionaries(self.param_dir.text())
            main = self.window()
            
            main.flags_dict = self.flags_dict
            main.params_dict = self.params_dict
            
            main.config_panel.flags_dict.clear()
            main.config_panel.flags_dict.update(self.flags_dict)
            
            main.config_panel.params_dict.clear()
            main.config_panel.params_dict.update(self.params_dict)
            self.flag_params = True
        except Exception as e:
            print(e)
            self.flag_params = False

    def toggle_rerun_fields(self, enabled: bool):
        self.config_file.setEnabled(enabled)
        self.config_file_button.setEnabled(enabled)
        self.dropdown.setEnabled(enabled)

    def select_config_file(self):
        path, _ = QFileDialog.getOpenFileName(self, dir=os.getcwd())
        if path:
            self.config_file.setText(path)
            self.load_config_file(path)

    def load_config_file(self, path):
        self.dropdown.clear()
        try:
            with open(path) as file:
                line = file.read().splitlines()[1]
                line = line.replace("order =", "").replace("_py", "").replace("_R", "")
                values = [el.strip() for el in line.split(",")]
                self.dropdown.addItems(values)
                self.dropdown.setEnabled(True)
        except Exception as e:
            print(e)

# ----------------------------
# Workflow list
# ----------------------------

class WorkflowListPanel(QWidget):
    def __init__(self):
        super().__init__()

        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.setSpacing(2)

        self.workflow_list = QListWidget()
        self.workflow_list.setDragDropMode(QListWidget.InternalMove)
        self.workflow_list.setSelectionMode(QAbstractItemView.SingleSelection)
        self.workflow_list.model().rowsMoved.connect(self.keep_general_on_top)

        general_step = WorkflowStep("General")
        general_item = QListWidgetItem("General")
        general_item.setData(Qt.UserRole, general_step)
        # make it non-removable / non-draggable
        general_item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled)
        
        self.workflow_list.addItem(general_item)

        self.remove_button = QPushButton("Remove Selected Step")
        self.save_config_button = QPushButton("Save configuration")
        self.load_config_button = QPushButton("Load configuration")
        self.write_disk_button = QPushButton("Generate Files")

        controls = QHBoxLayout()
        controls.addWidget(self.remove_button)
        
        bottom_buttons = QVBoxLayout()
        bottom_buttons.addWidget(self.save_config_button)
        bottom_buttons.addWidget(self.load_config_button)
        bottom_buttons.addWidget(self.write_disk_button)

        self.layout.addLayout(controls)
        self.layout.addWidget(self.workflow_list)
        self.layout.addLayout(bottom_buttons)

        self.block_list = None

        self.remove_button.clicked.connect(self.remove_step)
        self.save_config_button.clicked.connect(self.save_configuration)
        self.load_config_button.clicked.connect(self.load_configuration)
        self.write_disk_button.clicked.connect(self.generate_files)

        self.workflow_list.keyPressEvent = self.handle_key_press

    def set_block_source(self, block_list):
        self.block_list = block_list

    def add_step(self):
        if not self.block_list:
            return

        item = self.block_list.currentItem()
        if not item:
            return

        step = WorkflowStep(item.text())

        list_item = QListWidgetItem(step.name)
        list_item.setData(Qt.UserRole, step)

        self.workflow_list.addItem(list_item)
        self.workflow_list.setCurrentItem(list_item)

    def remove_step(self):
        item = self.workflow_list.currentItem()
    
        if item is None:
            return
    
        if item.text() == "General":
            return  # protect General
    
        row = self.workflow_list.row(item)
        self.workflow_list.takeItem(row)

    def handle_key_press(self, event):
        if event.key() == Qt.Key_Delete:
            self.remove_step()
        else:
            QListWidget.keyPressEvent(self.workflow_list, event)
    
    def keep_general_on_top(self, *args):
        for i in range(self.workflow_list.count()):
            item = self.workflow_list.item(i)
            if item.text() == "General" and i != 0:
                self.workflow_list.takeItem(i)
                self.workflow_list.insertItem(0, item)
                break
    
    def get_configuration_dict(self):
        workflow_data = {
            "global_settings": {
                "ms_file": self.main_window.global_settings.ms_file.text(),
                "meta_file": self.main_window.global_settings.meta_file.text(),
                "param_dir": self.main_window.global_settings.param_dir.text(),
                "work_dir": self.main_window.global_settings.work_dir.text(),
                "config_file": self.main_window.global_settings.config_file.text(),
                "rerun": self.main_window.global_settings.rerun_checkbox.isChecked(),
                "rerun_start": self.main_window.global_settings.dropdown.currentText(),
                "docker": self.main_window.global_settings.docker_checkbox.isChecked()
            },
            "workflow": []
        }
    
        for i in range(self.workflow_list.count()):
            item = self.workflow_list.item(i)
            step = item.data(Qt.UserRole)
    
            workflow_data["workflow"].append({
                "name": step.name,
                "params": step.params
            })
    
        return workflow_data
    
    def save_configuration(self):

        default_dir = self.main_window.global_settings.work_dir.text()

        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Save Workflow Configuration",
            default_dir,
            "JSON Files (*.json)"
        )
            
        if not filename:
            return
    
        workflow_data = self.get_configuration_dict()
    
        with open(filename, "w") as f:
            json.dump(workflow_data, f, indent=2)
    
    
    def load_configuration(self):
    
        filename, _ = QFileDialog.getOpenFileName(
            self,
            "Load Workflow Configuration",
            "",
            "JSON Files (*.json)"
        )
    
        if not filename:
            return
    
        with open(filename, "r") as f:
            workflow_data = json.load(f)
            
            settings = workflow_data.get("global_settings", {})
    
            gs = self.main_window.global_settings
            
            gs.ms_file.setText(settings.get("ms_file", ""))
            gs.meta_file.setText(settings.get("meta_file", ""))
            gs.param_dir.setText(settings.get("param_dir", ""))
            gs.work_dir.setText(settings.get("work_dir", ""))
            
            gs.rerun_checkbox.setChecked(settings.get("rerun", False))
            gs.config_file.setText(settings.get("config_file", ""))
            config_path = settings.get("config_file", "")
            if config_path:
                gs.load_config_file(config_path)
            gs.dropdown.setCurrentText(settings.get("rerun_start", ""))
            gs.docker_checkbox.setChecked(settings.get("docker", True))
            
            # automatically load data
            gs.validate()
    
        self.workflow_list.clear()
    
        for step_data in workflow_data["workflow"]:
    
            step = WorkflowStep(step_data["name"])
            step.params = step_data.get("params", {})
            item = QListWidgetItem(step.name)
            item.setData(Qt.UserRole, step)
    
            # protect General item
            if step.name == "General":
                item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled)
    
            self.workflow_list.addItem(item)
    
        self.keep_general_on_top()
        
    def generate_files(self):

        gs = self.main_window.global_settings
        workflow_dict = self.get_configuration_dict()
        workdir = os.path.abspath(gs.work_dir.text())
        data_fp = os.path.abspath(gs.ms_file.text())
        metadata_fp = os.path.abspath(gs.meta_file.text())
        rerun = gs.rerun_checkbox.isChecked()
        if rerun:
            rerun_start = gs.dropdown.currentText()
            rerun_counter = int(gs.dropdown.itemText(gs.dropdown.count() - 1).split("_")[0])
        else:
            rerun_start = None
            rerun_counter = 0
        flags_dict = gs.flags_dict
        params_dict = gs.params_dict
        docker_setup = gs.docker_checkbox.isChecked()
        try:
            gui_funcs.write_output(workflow_dict["workflow"], workdir, data_fp, metadata_fp, flags_dict, params_dict, rerun, rerun_start, rerun_counter, docker_setup)
            gs.status_label.setText("Files successfully generated ✅")

            QTimer.singleShot(
                5000,
                lambda: gs.status_label.clear()
            )
            
        except:
            gs.status_label.setText("Failed generating files ❌")

# ----------------------------
# Main Window
# ----------------------------
class MainWindow(QMainWindow):
    def __init__(self, flags_dict=None, params_dict=None, context=None):
        super().__init__()

        self.setWindowTitle("Parameter Generation")
        self.resize(1200, 700)

        self.context = context or {
            "metadata": {},
            "data_header": []
        }
        self.flags_dict = flags_dict or {}
        self.params_dict = params_dict or {}

        central = QWidget()
        grid = QGridLayout(central)
        grid.setSpacing(6)

        # Left top: Workflow generator
        self.block_list = QListWidget()
        self.block_list.addItems([
            "Filtering_long", "Filtering_wide", "Normalisation",
            "Imputation", "Correction", "Summarisation"
        ])

        self.add_button = QPushButton("Add Selected Step")

        self.left_top_box = QGroupBox("Available Modules")
        left_layout = QVBoxLayout(self.left_top_box)
        left_layout.setContentsMargins(5, 5, 5, 5)
        left_layout.setSpacing(4)

        left_layout.addWidget(self.block_list)
        left_layout.addWidget(self.add_button)

        # Right top: Global settings
        self.global_settings = GlobalSettingsPanel(self.context)
        self.global_settings.rerun_checkbox.stateChanged.connect(self.update_blocks)

        self.settings_box = QGroupBox("Settings")
        settings_layout = QVBoxLayout(self.settings_box)
        settings_layout.setContentsMargins(5, 5, 5, 5)

        settings_layout.addWidget(self.global_settings)

        #left bottom: Workflow
        self.workflow_panel = WorkflowListPanel()
        self.workflow_panel.main_window = self

        self.workflow_box = QGroupBox("Workflow")
        workflow_layout = QVBoxLayout(self.workflow_box)
        workflow_layout.setContentsMargins(5, 5, 5, 5)

        workflow_layout.addWidget(self.workflow_panel)

        # right bottom: Module configurator
        self.config_panel = ConfigPanel(self.flags_dict, self.params_dict, self.context)

        self.config_box = QGroupBox("Configuration")
        config_layout = QVBoxLayout(self.config_box)
        config_layout.setContentsMargins(5, 5, 5, 5)

        config_layout.addWidget(self.config_panel)

        grid.addWidget(self.left_top_box, 0, 0)
        grid.addWidget(self.settings_box, 0, 1)
        grid.addWidget(self.workflow_box, 1, 0)
        grid.addWidget(self.config_box, 1, 1)

        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 9)

        grid.setRowStretch(0, 1)
        grid.setRowStretch(1, 9)

        self.setCentralWidget(central)

        # Connections
        self.workflow_panel.set_block_source(self.block_list)

        self.add_button.clicked.connect(self.workflow_panel.add_step)

        self.workflow_panel.workflow_list.currentItemChanged.connect(
            self.config_panel.load_config_for_step
        )

        self.global_settings.settings_valid.connect(self.toggle_ui)

        # Initial state
        self.add_button.setEnabled(False)
        self.block_list.setEnabled(False)
        self.workflow_panel.workflow_list.setEnabled(False)
        self.workflow_panel.remove_button.setEnabled(False)
        self.workflow_panel.write_disk_button.setEnabled(False)
        self.workflow_panel.save_config_button.setEnabled(False)
        self.config_panel.setEnabled(False)

        # Enable Add button only if block selected
        self.block_list.currentItemChanged.connect(
            lambda: self.add_button.setEnabled(self.block_list.currentItem() is not None)
        )

        # Enable Remove button only if workflow item selected
        self.workflow_panel.remove_button.setEnabled(False)
        self.workflow_panel.workflow_list.currentItemChanged.connect(
            lambda current, _: self.workflow_panel.remove_button.setEnabled(current is not None)
        )
        
        general_item = self.workflow_panel.workflow_list.item(0)
        self.workflow_panel.workflow_list.setCurrentItem(general_item)
        
    # Enable UI after validation
    def toggle_ui(self, enabled):
        self.block_list.setEnabled(enabled)
        self.workflow_panel.workflow_list.setEnabled(enabled)
        self.workflow_panel.remove_button.setEnabled(enabled)
        self.workflow_panel.write_disk_button.setEnabled(enabled)
        self.workflow_panel.save_config_button.setEnabled(enabled)
        self.config_panel.setEnabled(enabled)

    # Update block list based on rerun
    def update_blocks(self, state):
        names = [self.block_list.item(i).text() for i in range(self.block_list.count())]

        if state:
            for i in range(self.block_list.count()):
                if self.block_list.item(i).text() == "Filtering_long":
                    self.block_list.takeItem(i)
                    break
        else:
            if "Filtering_long" not in names:
                self.block_list.addItem("Filtering_long")
        
# ----------------------------
# Run App
# ----------------------------
if __name__ == "__main__":
    # Use existing QApplication if present
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)

    window = MainWindow()
    window.show()

    if not QApplication.instance().startingUp():  
        sys.exit(app.exec())