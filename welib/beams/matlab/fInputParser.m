classdef fInputParser < handle
  properties
    %% FIXME: set input checking for these properties
    CaseSensitive = false;
    FunctionName  = '';
    KeepUnmatched = false;
    PartialMatching = false; % FIXME: unimplemented (and default should be true)
    StructExpand    = true;
  end

  properties (SetAccess = protected)
    Parameters    ;%= cell ();
    Results       = struct ();
    Unmatched     = struct ();
    UsingDefaults ;%= cell ();
  end

  properties (Access = protected)
    %% Since Required and Optional are ordered, they get a cell array of
    %% structs with the fields 'name', 'def' (default), and 'val' (validator).
    Required ;%= cell ();
    Optional ;%= cell ();
    %% Parameter and Switch are unordered so we have a struct whose fieldnames
    %% are the argname, and values are a struct with fields 'def' and 'val'
    Parameter = struct ();
    Switch    = struct ();

    %% List of Parameter and Switch names to ease searches
    ParameterNames ;%= cell ();
    SwitchNames    ;%= cell ();

    %% When checking for fieldnames in a Case Insensitive way, this variable
    %% holds the correct identifier for the last searched named using the
    %% is_argname method.
    last_name = '';
  end

  properties (Access = protected, Constant = true)
    %% Default validator, always returns scalar true.
    def_val = @(x) true;
  end

  methods
    function set.PartialMatching (this, val)
      if (val)
        error ('fInputParser: PartialMatching is not yet implemented');
      end
    end

    function addRequired (this, name, val)
      if ~exist('val','var'); val=fInputParser.def_val; end;

      if (nargin < 2 || nargin > 3)
        print_usage ();
      elseif (numel (this.Optional) || numfields (this.Parameter) || numfields (this.Switch))
        error ('fInputParser.addRequired: can''t have a Required argument after Optional, Parameter, or Switch');
      end
      this.validate_name ('Required', name);
      this.Required{end+1} = struct ('name', name, 'val', val);
    end

    function addOptional (this, name, def, val)
      if ~exist('val','var'); val=fInputParser.def_val; end;

      if (nargin < 3 || nargin > 4)
        print_usage ();
      elseif (numfields (this.Parameter) || numfields (this.Switch))
        error ('fInputParser.Optional: can''t have Optional arguments after Parameter or Switch');
      end
      this.validate_name ('Optional', name);
      this.Optional{end+1} = struct ('name', name, 'def', def, 'val', val);
    end

    function addParamValue (this, name, def, val)
      if ~exist('val','var'); val=fInputParser.def_val; end;

      if (nargin < 3 || nargin > 4)
        print_usage ();
      end
      this.addParameter (name, def, val);
    end

    function addParameter (this, name, def, varargin)
      if (nargin < 3 || nargin > 6)
        print_usage ();
      end

      n_opt = numel (varargin);

      if (n_opt == 0 || n_opt == 2)
        val = fInputParser.def_val;
      else % n_opt is 1 or 3
        val = varargin{1};
      end

      if (n_opt == 0 || n_opt == 1)
          match_priority = 1;
      else % n_opt is 2 or 3
          if (~ strcmpi (varargin{end-1}, 'PartialMatchPriority'))
              error ('fInputParser.addParameter: unrecognized option');
          end
          match_priority = varargin{end};
          validateattributes (match_priority, {'numeric'}, {'positive', 'integer'}, 'fInputParser.addParameter', 'PartialMatchPriority');
      end

      this.validate_name ('Parameter', name);
      this.Parameter.(name).def = def;
      this.Parameter.(name).val = val;
  end

    function addSwitch (this, name)
      if (nargin ~= 2)
        print_usage ();
      end
      this.validate_name ('Switch', name);
      this.Switch.(name).def = false;
    end

    function parse (this, varargin)
      this.Results       = struct();
      this.Unmatched     = struct();
      this.UsingDefaults = cell (0);
      if (numel (varargin) < numel (this.Required))
        if (this.FunctionName)
          print_usage (this.FunctionName);
        else
          this.error ('fInputParser.parse: not enough input arguments');
        end
      end
      pnargin = numel(varargin);

      this.ParameterNames = fieldnames (this.Parameter);
      this.SwitchNames    = fieldnames (this.Switch);

      %% Evaluate the Required arguments first
      nReq = numel (this.Required);
      for idx = 1:nReq
        req = this.Required{idx};
        this.validate_arg (req.name, req.val, varargin{idx});
      end

      vidx = nReq;  % current index in varargin

      %% Search for a list of Optional arguments
      idx  = 0;     % current index on the array of Optional
      nOpt = numel(this.Optional);
      while (vidx < pnargin && idx < nOpt)
        opt = this.Optional{idx+1};
        in  = varargin{vidx+1};
        if ((this.is_argname ('Parameter', in) && vidx < pnargin) || this.is_argname ('Switch', in))
          %% This looks like an optional parameter/value pair or a
          %% switch, not an positional option.  This does mean that
          %% positional options cannot be strings named like parameter
          %% keys.  See bug %50752.
          idx = idx-1;
          vidx = vidx-1;
          break
        end
        try
          valid_option = opt.val (in);
        catch
          valid_option = false;
        end
        if (~ valid_option)
          %% If it does not match there's two options:
          %%    1) input is actually wrong and we should error;
          %%    2) it's a Parameter or Switch name and we should use
          %%       the default for the rest;
          %%    3) it's a struct with the Parameter pairs.
          if (ischar (in) || (this.StructExpand && isstruct (in) && isscalar (in)))
            idx = idx-1;
            vidx = vidx-1;
            break
          else
            this.error (sprintf (['failed validation of %s\n',  'Validation function: %s'],  (opt.name), disp(opt.val)));
          end
        end
        this.Results.(opt.name) = in;
      end

      %% Fill in with defaults of missing Optional
      while (idx < nOpt)
        idx=idx+1;
        opt = this.Optional{idx};
        this.UsingDefaults{end+1} = opt.name;
        this.Results.(opt.name) = opt.def;
      end

      %% Search unordered Options (Switch and Parameter)
      while (vidx < pnargin)
        vidx=vidx+1;
        name = varargin{vidx};

        if (this.StructExpand && isstruct (name) && isscalar (name))
          expanded_options = [fieldnames(name) struct2cell(name)]';
          expanded_options = expanded_options(:);
          n_new_args = numel (expanded_options) -1;
          pnargin = pnargin + n_new_args;
          varargin(vidx+n_new_args+1:pnargin) = varargin(vidx+1:end);
          varargin(vidx:vidx+n_new_args) = expanded_options;
          name = varargin{vidx};
        end

        if (~ ischar (name))
          this.error ('non-string for Parameter name or Switch');
        end

        if (this.is_argname ('Parameter', name))
          if (vidx+1 > pnargin)
            this.error (sprintf ('no matching value for option %s', name));
          end
          %keyboard
          this.validate_arg (this.last_name,this.Parameter.(this.last_name).val, varargin{vidx+1});
        elseif (this.is_argname ('Switch', name))
          this.Results.(this.last_name) = true;
        else
          if (vidx+1 <= pnargin && this.KeepUnmatched)
            this.Unmatched.(name) = varargin{vidx+1};
          else
            this.error (sprintf ('argument %s is not a valid parameter',  name));
          end
        end
        vidx=vidx+1;
      end % while
      %% Add them to the UsingDefaults list
      this.add_missing ('Parameter');
      this.add_missing ('Switch');
      %%

    end

    function disp (this)
      if (nargin ~= 1)
        print_usage ();
      end
      fprintf ('fInputParser object with properties:\n\n');
      b2s = @(x) any (x);
      fprintf (['   CaseSensitive   : %s\n   FunctionName    : %s\n' ...
               '   KeepUnmatched   : %s\n   PartialMatching : %s\n' ...
               '   StructExpand    : %s\n\n'], b2s (this.CaseSensitive), b2s (this.FunctionName), b2s (this.KeepUnmatched), b2s (this.PartialMatching), b2s (this.StructExpand));
      if isempty(this.Parameters) 
          fprintf ('Defined parameters:\n\n   {%s}\n', '');
      else
          fprintf ('Defined parameters:\n\n   {%s}\n', strjoin (this.Parameters, ', '));
      end
    end
  end

  methods (Access = private)
    function validate_name (this, type, name)
      if (~ isvarname (name))
        error ('fInputParser.add%s: NAME is an invalid identifier', method);
      elseif (any (strcmpi (this.Parameters, name)))
        %% Even if CaseSensitive is 'on', we still shouldn't allow
        %% two args with the same name.
        error ('fInputParser.add%s: argname %s has already been specified', type, name);
      end
      this.Parameters{end+1} = name;
    end

    function validate_arg (this, name, val, in)
        if (~val(in))
%             keyboard
          this.error (sprintf ('failed validation of %s with %s',  name, func2str (val)));
        end
        this.Results.(name) = in;
    end

    function r = is_argname (this, type, name)
      if (this.CaseSensitive)
        r = isfield (this.(type), name);
        this.last_name = name;
      else
        fnames = this.([type 'Names']);
        l = strcmpi (name, fnames);
        r = any (l(:));
        if (r)
          this.last_name = fnames{l};
        end
      end
    end

    function add_missing (this, type)
      unmatched = setdiff (fieldnames (this.(type)), fieldnames (this.Results));
      for namec = unmatched(:)'
        name = namec{1};
        this.UsingDefaults{end+1} = name;
        this.Results.(name) = this.(type).(name).def;
      end
    end

    function error (this, msg)
      where = '';
      if (this.FunctionName)
        where = [this.FunctionName ': '];
      end
      error ('%s%s', where, msg);
    end
  end

end

%!function p = create_p ()
%!  p = fInputParser ();
%!  p.CaseSensitive = true;
%!  p.addRequired ('req1', @(x) ischar (x));
%!  p.addOptional ('op1', 'val', @(x) any (strcmp (x, {'val', 'foo'})));
%!  p.addOptional ('op2', 78, @(x) x > 50);
%!  p.addSwitch ('verbose');
%!  p.addParameter ('line', 'tree', @(x) any (strcmp (x, {'tree', 'circle'})));
%!end

%% check normal use, only required are given
%!test
%! p = create_p ();
%! p.parse ('file');
%! r = p.Results;
%! assert (r.req1, 'file');
%! assert (sort (p.UsingDefaults), sort ({'op1', 'op2', 'verbose', 'line'}));
%! assert ({r.req1, r.op1, r.op2, r.verbose, r.line},
%!         {'file', 'val', 78,    false,     'tree'});

%% check normal use, but give values different than defaults
%!test
%! p = create_p ();
%! p.parse ('file', 'foo', 80, 'line', 'circle', 'verbose');
%! r = p.Results;
%! assert ({r.req1, r.op1, r.op2, r.verbose, r.line},
%!         {'file', 'foo', 80,    true,      'circle'});

%% check optional is skipped and considered Parameter if unvalidated string
%!test
%! p = create_p ();
%! p.parse ('file', 'line', 'circle');
%! r = p.Results;
%! assert ({r.req1, r.op1, r.op2, r.verbose, r.line},
%!         {'file', 'val', 78,    false,     'circle'});

%% check case insensitivity
%!test
%! p = create_p ();
%!  p.CaseSensitive = false;
%! p.parse ('file', 'foo', 80, 'LiNE', 'circle', 'vERbOSe');
%! r = p.Results;
%! assert ({r.req1, r.op1, r.op2, r.verbose, r.line},
%!         {'file', 'foo', 80,    true,      'circle'});

%% check KeepUnmatched
%!test
%! p = create_p ();
%! p.KeepUnmatched = true;
%! p.parse ('file', 'foo', 80, 'line', 'circle', 'verbose', 'extra', 50);
%! assert (p.Unmatched.extra, 50);

%% check error when missing required
%!error <not enough input arguments>
%! p = create_p ();
%! p.parse ();

%% check error when given required does not validate
%!error <failed validation of >
%! p = create_p ();
%! p.parse (50);

%% check error when given optional does not validate
%!error <is not a valid parameter>
%! p = create_p ();
%! p.parse ('file', 'no-val');

%% check error when given Parameter does not validate
%!error <failed validation of >
%! p = create_p ();
%! p.parse ('file', 'foo', 51, 'line', 'round');

%% check alternative method (obj, ...) API
%!function p2 = create_p2 ();
%!  p2 = fInputParser;
%!  addRequired (p2, 'req1', @(x) ischar (x));
%!  addOptional (p2, 'op1', 'val', @(x) any (strcmp (x, {'val', 'foo'})));
%!  addOptional (p2, 'op2', 78, @(x) x > 50);
%!  addSwitch (p2, 'verbose');
%!  addParameter (p2, 'line', 'tree', @(x) any (strcmp (x, {'tree', 'circle'})));
%!end

%% check normal use, only required are given
%!test
%! p2 = create_p2 ();
%! parse (p2, 'file');
%! r = p2.Results;
%! assert ({r.req1, r.op1, r.op2, r.verbose, r.line},
%!         {'file', 'val', 78,    false,     'tree'});
%! assert (sort (p2.UsingDefaults), sort ({'op1', 'op2', 'verbose', 'line'}));

%% check normal use, but give values different than defaults
%!test
%! p2 = create_p2 ();
%! parse (p2, 'file', 'foo', 80, 'line', 'circle', 'verbose');
%! r = p2.Results;
%! assert ({r.req1, r.op1, r.op2, r.verbose, r.line},
%!         {'file', 'foo', 80,    true,      'circle'});

%% We must not perform validation of default values
%!test <*45837>
%! p = fInputParser;
%! p.addParameter ('Dir', [], @ischar);
%! p.parse ();
%! assert (p.Results.Dir, []);

%!test
%! p = fInputParser;
%! p.addParameter ('positive', -1, @(x) x > 5);
%! p.parse ();
%! assert (p.Results.positive, -1);

%% Throw an error on validation of optional argument to check that it
%% is caught without preventing continuation into param/value pairs.
%!test
%! p = fInputParser ();
%! p.addOptional ('err', 'foo', @error);
%! p.addParameter ('not_err', 'bar', @ischar);
%! p.parse ('not_err', 'qux');
%! assert (p.Results.err, 'foo')
%! assert (p.Results.not_err, 'qux')


%% With more Parameters to test StructExpand
%!function p3 = create_p3 ();
%!  p3 = fInputParser;
%!  addOptional (p3, 'op1', 'val', @(x) any (strcmp (x, {'val', 'foo'})));
%!  addOptional (p3, 'op2', 78, @(x) x > 50);
%!  addSwitch (p3, 'verbose');
%!  addParameter (p3, 'line', 'tree', @(x) any (strcmp (x, {'tree', 'circle'})));
%!  addParameter (p3, 'color', 'red', @(x) any (strcmp (x, {'red', 'green'})));
%!  addParameter (p3, 'style', 'tt', @(x) any (strcmp (x, {'tt', 'f', 'i'})));
%!end

%% Test StructExpand
%!test
%! p3 = create_p3 ();
%! p3.parse (struct ('line', 'circle', 'color', 'green'));
%! assert (p3.Results, struct ('op1', 'val', 'op2', 78, 'verbose', false,
%!                             'line', 'circle', 'color', 'green',
%!                             'style', 'tt'))

%!test
%! p3 = create_p3 ();
%! p3.parse (struct ('line', 'circle', 'color', 'green'), 'line', 'tree');
%! assert (p3.Results.line, 'tree')
%! p3.parse ('line', 'tree', struct ('line', 'circle', 'color', 'green'));
%! assert (p3.Results.line, 'circle')

%!test % unmatched parameters with StructExpand
%! p3 = create_p3 ();
%! p3.KeepUnmatched = true;
%! p3.parse (struct ('line', 'circle', 'color', 'green', 'bar', 'baz'));
%! assert (p3.Unmatched.bar, 'baz')

%% The validation for the second optional argument throws an error with
%% a struct so check that we can handle it.
%!test
%! p3 = create_p3 ();
%! p3.parse ('foo', struct ('color', 'green'), 'line', 'tree');
%! assert (p3.Results.op1, 'foo')
%! assert (p3.Results.line, 'tree')
%! assert (p3.Results.color, 'green')
%! assert (p3.Results.verbose, false)


%% Some simple tests for addParamValue since all the other ones use add
%% addParameter but they use the same codepath.
%!test
%! p = fInputParser;
%! addParameter (p, 'line', 'tree', @(x) any (strcmp (x, {'tree', 'circle'})));
%! addParameter (p, 'color', 'red', @(x) any (strcmp (x, {'red', 'green'})));
%! p.parse ('line', 'circle');
%! assert ({p.Results.line, p.Results.color}, {'circle', 'red'})

%!test
%! p = fInputParser;
%! p.addParameter ('foo', 'bar', @ischar);
%! p.parse ();
%! assert (p.Results, struct ('foo', 'bar'))
%! p.parse ('foo', 'qux');
%! assert (p.Results, struct ('foo', 'qux'))

%% This behaviour means that a positional option can never be a string
%% that is the name of a parameter key.  This is required for Matlab
%% compatibility.
%!test <*50752>
%! p = fInputParser ();
%! p.addOptional ('op1', 'val');
%! p.addParameter ('line', 'tree');
%! p.parse ('line', 'circle');
%! assert (p.Results, struct ('op1', 'val', 'line', 'circle'))
%!
%! p = fInputParser ();
%! p.addOptional ('op1', 'val1');
%! p.addOptional ('op2', 'val2');
%! p.addParameter ('line', 'tree');
%! p.parse ('line', 'circle');
%! assert (p.Results.op1, 'val1')
%! assert (p.Results.op2, 'val2')
%! assert (p.Results.line, 'circle')
%!
%! %% If there's enough arguments to fill the positional options and
%! %% param/key, it still skips positional options.
%! p = fInputParser ();
%! p.addOptional ('op1', 'val1');
%! p.addOptional ('op2', 'val2');
%! p.addParameter ('line', 'tree');
%! p.parse ('line', 'circle', 'line', 'rectangle');
%! assert (p.Results, struct ('op1', 'val1', 'op2', 'val2',
%!                            'line', 'rectangle'))
%!
%! %% Even if the key/param fails validation, it does not backtrack to
%! %% check if the values are valid positional options.
%! p = fInputParser ();
%! p.addOptional ('op1', 'val1', @ischar);
%! p.addOptional ('op2', 'val2', @isnumeric);
%! p.addParameter ('line', 'circle', @ischar);
%! fail ('p.parse ('line', 89)', 'failed validation of LINE')
%!
%! p = fInputParser ();
%! p.addOptional ('op1', 'val1');
%! p.addParamValue ('line', 'circle', @ischar);
%! fail ('p.parse ('line', 'line', 89)',
%!       'non-string for Parameter name or Switch')

%!test <*50752>
%! %% This fails in Matlab but works in Octave.  It is a bug there
%! %% that we do not replicate.
%! p = fInputParser ();
%! p.addOptional ('op1', 'val1');
%! p.addParameter ('line', 'circle');
%! p.parse ('line');
%! assert (p.Results, struct ('op1', 'line', 'line', 'circle'))

%!test <*50752>
%! p = fInputParser;
%! p.addOptional ('op1', 'val1');
%! p.addSwitch ('line');
%! p.parse ('line');
%! assert (p.Results.op1, 'val1')
%! assert (p.Results.line, true)

%!test
%! p = fInputParser;
%! p.addParameter ('a', []);
%! p.addParameter ('b', []);
%! p.parse ('a', 1);
%! p.parse ('b', 1);
%! assert (p.Results, struct ('a', [], 'b', 1));
%! assert (p.UsingDefaults, {'a'});

%!test
%! p = fInputParser;
%! p.addParameter ('b', []);
%! p.KeepUnmatched = true;
%! p.parse ('a', 1);
%! p.parse ('b', 1);
%! assert (p.Results, struct ('b', 1));
%! assert (p.Unmatched, struct ());

%% Test for patch %9241
%!error<failed validation of A with ischar>
%! p = fInputParser;
%! p.addParameter ('a', [], @ischar);
%! p.parse ('a', 1);
