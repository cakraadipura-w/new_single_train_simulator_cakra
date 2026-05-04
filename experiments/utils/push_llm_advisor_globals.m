function push_llm_advisor_globals(enable_flag, profile)
%PUSH_LLM_ADVISOR_GLOBALS Set LLM advisor globals on each MATLAB session.
%
% Usage:
%   push_llm_advisor_globals(true, profile_struct)
%   wait(parfevalOnAll(@push_llm_advisor_globals, 0, true, profile_struct));

    global llm_advisor_enabled llm_advisor_profile

    if nargin < 1 || isempty(enable_flag)
        llm_advisor_enabled = false;
    else
        llm_advisor_enabled = logical(enable_flag);
    end

    if nargin < 2 || isempty(profile)
        llm_advisor_profile = struct();
    else
        llm_advisor_profile = profile;
    end
end