import sys
import traceback 


try:
    # see https://github.com/nir0s/backtrace
    import backtrace
    from backtrace import Fore, Style
    STYLES = {
        'backtrace': Fore.YELLOW + '{0}',
        'error': Fore.RED + Style.BRIGHT + '{0}',
        'line': Fore.RED + Style.BRIGHT + '{0:4s}',
#         'module': '{0}',
        'module': '{0:20s}',
#         'context': Style.BRIGHT + Fore.GREEN + '{0}',
        'context': '',
        'call': Fore.YELLOW + '> ' + Style.BRIGHT + '{0}',
    }


    backtrace.hook(
        reverse=False,
        align=False,
        strip_path=True,
        enable_on_envvar_only=False,
        on_tty=False,
        conservative=False,
        styles=STYLES)
except:
    print('Cannot use clean_exceptions, module `backtrace` not installed')

    # more code...

#     # if you wanna restore the default hook...
#     backtrace.unhook()


# def MyExceptionHook(etype, value, trace):
#     """
#     Handler for all unhandled exceptions.
#     :param `etype`: the exception type (`SyntaxError`, `ZeroDivisionError`, etc...);
#     :type `etype`: `Exception`
#     :param string `value`: the exception error message;
#     :param string `trace`: the traceback header, if any (otherwise, it prints the
#      standard Python header: ``Traceback (most recent call last)``.
# 
# 
#     Documentaion of traceback.print_exceptions 
#         This differs from print_tb() in the following ways: (1) if
#         traceback is not None, it prints a header "Traceback (most recent
#         call last):"; (2) it prints the exception type and value after the
#         stack trace; (3) if type is SyntaxError and value has the
#         appropriate format, it prints the line where the syntax error
#         occurred with a caret on the next line indicating the approximate
#         position of the error.
# 
#     """
#     outputstream = sys.stderr
#     # Printing exception
#     print(traceback.__file__)
#     print('>>> Trace:')
#     #traceback.print_exception(etype, value, trace)
# 
#     
# 
#     for line in traceback.TracebackException(type(value), value, trace, limit=None).format(chain=True):
# #         print(line, file=sys.stderr, end="")
#         print(type(line))
# 
# 
# 
#     # Then showing to user the last error
#     #frame = wx.GetApp().GetTopWindow()
#     tmp = traceback.format_exception(etype, value, trace)
#     print('>>> Exception:')
#     print('The following exception occured:\n\n'+ tmp[-1]  + '\n'+tmp[-2].strip())
# 
# sys.excepthook = MyExceptionHook
