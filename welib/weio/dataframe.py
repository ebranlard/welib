import pandas as pd
from .converters import writeDataFrameToFormat, writeDataFrameAutoFormat


class WEIODataFrame(pd.DataFrame):
    """
    Custom DataFrame subclass providing:
     - format export capabilities
     - case insensitive columns
    and flexible initialization from a DataFrame or raw data/columns.
    """

    def __init__(self, data=None, index=None, columns=None, dtype=None, copy=None, cols_readonly=False):
        if isinstance(data, pd.DataFrame):
            # If initialized directly with a DataFrame, copy its data & index
            super().__init__(
                data=data.values,
                index=index if index is not None else data.index,
                columns=columns if columns is not None else data.columns,
                dtype=dtype,
                copy=copy,
            )
        else:
            super().__init__(
                data=data, index=index, columns=columns, dtype=dtype, copy=copy
            )
        self._cols_readonly = cols_readonly

    @property
    def _constructor(self):
        """Ensures pandas operations return a WEIODataFrame instance."""
        return WEIODataFrame

    @property
    def cols_readonly(self):
        return getattr(self, "_cols_readonly", False)

    @cols_readonly.setter
    def cols_readonly(self, value):
        self._cols_readonly = bool(value)

    def _get_case_insensitive_key(self, key):
        """Helper to resolve a string key case-insensitively against columns."""
        if isinstance(key, str) and key not in self.columns:
            key_lower = key.lower()
            mapping = {
                col.lower(): col for col in self.columns if isinstance(col, str)
            }
            if key_lower in mapping:
                return mapping[key_lower]
        return key

    def __setitem__(self, key, value):
        if not self.cols_readonly:
            super().__setitem__(key, value)
        else:
            # Resolving case-insensitive match for existing columns
            resolved_key = self._get_case_insensitive_key(key) if isinstance(key, str) else key
            
            # Allow modification of existing columns, but block adding new ones
            if isinstance(resolved_key, str) and resolved_key not in self.columns:
                raise TypeError(f"Cannot add new column '{key}'. Columns are read-only.")
            elif isinstance(key, list):
                new_cols = [k for k in key if (self._get_case_insensitive_key(k) if isinstance(k, str) else k) not in self.columns]
                if new_cols:
                    raise TypeError(f"Cannot add new columns {new_cols}. Columns are read-only.")
            
            super().__setitem__(resolved_key, value)

    def __getitem__(self, key):
        """Override bracket indexing to support case-insensitive column lookups."""
        if isinstance(key, str):
            resolved_key = self._get_case_insensitive_key(key)
            return super().__getitem__(resolved_key)
        elif isinstance(key, list):
            resolved_keys = [
                self._get_case_insensitive_key(k) if isinstance(k, str) else k
                for k in key
            ]
            return super().__getitem__(resolved_keys)
        return super().__getitem__(key)

    def __delitem__(self, key):
        if self.cols_readonly:
            raise TypeError(f"Cannot delete column '{key}'. Columns are read-only.")
        resolved_key = self._get_case_insensitive_key(key) if isinstance(key, str) else key
        super().__delitem__(resolved_key)


    def to_format(self, filename, fformat):
        """Write DataFrame to disk using specified format."""
        writeDataFrameToFormat(self, filename, fformat)

    def export(self, filename, fformat=None):
        """Write DataFrame to disk using specified format based on extension.

         .csv :   CSV
         .outb :  FASTOutputFile
         .pq, parquet :  Parquet

        """
        writeDataFrameAutoFormat(self, filename, fformat)

    def to_outb(self, filename):
        """Write DataFrame to FAST binary output (.outb) format."""
        writeDataFrameToFormat(self, filename, "outb")

    def to_parquet(self, filename, **kwargs):
        """Write DataFrame to Parquet format."""
        super().to_parquet(path=filename, **kwargs)
        #writeDataFrameToFormat(self, filename, "parquet") # will trigger a recursion

    def to_csv(self, filename, index=False, sep=',', **kwargs):
        """Write DataFrame to CSV using writeDataFrameToFormat handler."""
        super().to_csv(filename, sep=sep, index=index, **kwargs)
        #writeDataFrameToFormat(self, filename, "csv") # will trigger a recursion


