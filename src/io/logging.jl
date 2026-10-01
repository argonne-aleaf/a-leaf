# Copyright © 2026, UChicago Argonne, LLC. All Rights Reserved.
# SPDX-License-Identifier: BSD-3-Clause  (see LICENSE)

# ALEAF logging macros: wrap the stdlib log macros and flush stdout/stderr.

using Logging

macro aleaf_info(args...)
    log_call = Expr(:macrocall, Symbol("@info"), __source__, esc.(args)...)
    return quote
        $(log_call)
        flush(stderr)
        flush(stdout)
    end
end

macro aleaf_warn(args...)
    log_call = Expr(:macrocall, Symbol("@warn"), __source__, esc.(args)...)
    return quote
        $(log_call)
        flush(stderr)
        flush(stdout)
    end
end

macro aleaf_error(args...)
    log_call = Expr(:macrocall, Symbol("@error"), __source__, esc.(args)...)
    return quote
        $(log_call)
        flush(stderr)
        flush(stdout)
        msg = join(string.($(esc.(args)...)), " ")
        error("ERROR: " * msg)
    end
end

macro aleaf_debug(args...)
    log_call = Expr(:macrocall, Symbol("@debug"), __source__, esc.(args)...)
    return quote
        $(log_call)
        flush(stderr)
        flush(stdout)
    end
end
