# command_line.jl
#
# Two small helpers for reading options of the form  --name=value  and switches of the form  --name
# from the command line of a script.
#
# Example:   julia --project=.. some_script.jl --tend=500 --quick
#            option("tend", 2000.0)  returns 500.0   (parsed with the type of the default value)
#            option("seed", 1)       returns 1       (option not given: the default is returned)
#            flag("quick")           returns true

"""
    option(name, default)

Value of the command-line option `--name=value`, converted to the type of `default`.
Returns `default` if the option is not given.
"""
function option(name::AbstractString, default)
    prefix = "--" * name * "="
    for argument in ARGS
        if startswith(argument, prefix)
            text = argument[length(prefix)+1:end]
            return convert_option(text, default)
        end
    end
    return default
end

convert_option(text, default::AbstractString) = String(text)
convert_option(text, default::Number) = parse(typeof(default), text)

"""
    flag(name)

`true` if the switch `--name` is given on the command line.
"""
flag(name::AbstractString) = ("--" * name) in ARGS
