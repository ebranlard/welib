import pandas as pd
from .converters import writeDataFrameToFormat, writeDataFrameAutoFormat


class WEIODataFrame(pd.DataFrame):
    """
    Custom DataFrame subclass providing:
     - format export capabilities
     - case insensitive columns
    and flexible initialization from a DataFrame or raw data/columns.
    """

    def __init__(self, data=None, index=None, columns=None, dtype=None, copy=None):
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

    @property
    def _constructor(self):
        """Ensures pandas operations return a WEIODataFrame instance."""
        return WEIODataFrame

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
        self.to_parquet(path=filename, **kwargs)
        #writeDataFrameToFormat(self, filename, "parquet")

    def to_csv(self, filename):
        """Write DataFrame to CSV using writeDataFrameToFormat handler."""
        self.to_csv(filename, sep=sep, index=index, **kwargs)
        writeDataFrameToFormat(self, filename, "csv")


