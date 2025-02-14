function [Diary_k, Diary_r, Diary_w, Diary_WeWu, Diary_aggregate_savings, Diary_distribution, Diary_distribution_g_e, Diary_distribution_g_u] = Diary()
    global N;
    Diary_k = zeros(N,1);
    Diary_r = zeros(N,1);
    Diary_w = cell(N, 1);
    Diary_distribution = cell(N,1);
    Diary_distribution_g_e = cell(N,1);
    Diary_distribution_g_u = cell(N,1);
    Diary_Job_Value = cell(N,1);
    Diary_Vacancy_Value = zeros(N,1);
    Diary_Firm_Profit = zeros(N,1);
    Diary_w_unadjusted = cell(N,1);
    Diary_Theta = zeros(N,1);
    Diary_kgap = zeros(N,1);
    Diary_lambda_u = zeros(N,1);
end
