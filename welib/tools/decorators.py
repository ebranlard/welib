from functools import wraps

def require_attrs(attrs, message=""):
    """
    Decorator to verify instance attributes are set before running a method.
    
    :param attrs: A list or tuple of attribute names to check.
    :param message: Optional additional context for the exception.

    usage:

        @require_attrs(['HD','fnd'], 'Optional message')
        def mymethod(self):
            pass
    """
    def decorator(func):
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            for attr in attrs:
                if getattr(self, attr, None) is None:
                    err_msg = f"Attribute '{attr}' must be set before calling '{func.__name__}'."
                    if message:
                        err_msg += f" {message}"
                    raise ValueError(err_msg)
            return func(self, *args, **kwargs)
        return wrapper
    return decorator
