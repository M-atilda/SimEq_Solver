defmodule SimEq.Gauss do
  import Matrix

  #@fn      [float] solve_gaussian_eliminate(%Matrix, [float], (bool))
  #@param  coefficient matrix
  #@param  inhomogenious vector
  #@param  use pivoting or not(option)

  @moduledoc """
  provides the direct method "gaussian"
    for solving simultenious equations

  pivot is choiced in partial pivoting
  complete pivoting is not yet implemeted

  for the complete pivoting
  the history of swapping rows
      for use suitable value as pivot
    is memorized as taples' list
  ex) swap line 1 and 3 => [{1, 3}]

  whether the simultinious equations is solvable
    is not treated at this module
  it assumes pre-process has already been done
    (its rank is same as the number of variables)
  """




  # helpers for forward elimination
  defp get_pivot_index(matrix, i, false) do
    #partial pivoting
    candidates = Enum.drop(get_row(matrix, i), i-1)
    pivot_num = max(Enum.max(candidates), abs(Enum.min(candidates)))
    i + Enum.find_index(candidates, &(abs(&1) == pivot_num))
  end
  defp get_pivot_index(%Matrix{line: l_m} = matrix, i, true) do
    #complete pivoting
    targets = for k <- i..l_m do
      Enum.drop(get_line(matrix, k), i-1) end
      |> Enum.map(fn l -> Enum.take(l, l_m-(i-1)) end)
    max_vector = Enum.map(targets,
      fn l -> max(Enum.max(l), abs(Enum.min(l))) end)
    pivot_num = Enum.max(max_vector)
    max_i = Enum.find_index(max_vector, &(&1 == pivot_num))
    max_j = Enum.find_index(Enum.at(targets, max_i), &(&1 == pivot_num))
    {i+max_i, i+max_j}
  end


  defp div_sub_lines([], _, _, matrix), do: matrix
  defp div_sub_lines([h|tail], j, i, matrix) do
    div_sub_lines(tail, j+1, i, mul_R_l(matrix, j, i, -h))
  end

  defp fe_helper(%Matrix{line: l_m} = matrix, i,
    _, _, acm) when l_m == i, do: {matrix, acm}
  defp fe_helper(%Matrix{line: l_m} = matrix, i,
    false, _, _) do
    pivot_num = get_factor(matrix, i, i)
    div_nums = for j <- (i+1)..l_m, do: get_factor(matrix, j, i) / pivot_num
    fe_helper(div_sub_lines(div_nums, i+1, i, matrix),
      i+1,
      false, false, [])
  end
  defp fe_helper(%Matrix{line: l_m} = matrix, i,
    true, is_comp, acm) do
    if is_comp do
      {pivot_i, pivot_j} = get_pivot_index(matrix, i, true)
      swaped_matrix = mul_P_l(matrix, i, pivot_i)
                      |> mul_P_r(i, pivot_j)
      pivot_num = get_factor(swaped_matrix, i, i)
      div_nums = for j <- (i+1)..l_m do
        get_factor(swaped_matrix, j, i) / pivot_num
      end
      fe_helper(div_sub_lines(div_nums, i+1, i, swaped_matrix),
                i+1,
                true, true,
                [{i, pivot_j}|acm])
    else
      pivot_id = get_pivot_index(matrix, i, false)
      swaped_matrix = mul_P_l(matrix, i, pivot_id)
      pivot_num = get_factor(swaped_matrix, i, i)
      div_nums = for j <- (i+1)..l_m do
        get_factor(swaped_matrix, j, i) / pivot_num
      end
      fe_helper(div_sub_lines(div_nums, i+1, i, swaped_matrix),
                i+1,
                true, false,
                acm)
    end
  end


  def forward_eliminate(%Matrix{line: l_m, row: r_n}, _, _)
    when l_m != r_n, do: raise "math error <forward_eliminate@gauss.ex>"
  def forward_eliminate(m, pivoting, is_comp) do
    fe_helper(m, 1, pivoting, is_comp, [])
  end



  # helpers for backward substitution
  defp calc_result([], result), do: result
  defp calc_result(c, []) do
    last_l = List.last c
    calc_result(Enum.slice(c, 0..-2), [List.last(last_l) / Enum.at(last_l, -2)])
  end
  defp calc_result(c, result) do
    [inhom_val|last_l] = Enum.reverse(List.last c)
    sub_val = Enum.zip(Enum.reverse(result), last_l)
              |> Enum.reduce(0, fn ({f, s}, acm) -> (f*s)+acm end)
    u_ii = Enum.at(last_l, length(result))
    calc_result(Enum.slice(c, 0..-2),
      [(inhom_val-sub_val) / u_ii|result])
  end

  defp modify_results_order(result, []), do: result
  defp modify_results_order(result, [{f_id, s_id}|tail]) do
    f_val = Enum.at(result, f_id-1)
    s_val = Enum.at(result, s_id-1)
    new_result = List.replace_at(result, f_id-1, s_val)
                 |> List.replace_at(s_id-1, f_val)
    modify_results_order(new_result, tail)
  end


  def backward_substitute(%Matrix{contents: c}, hist) do
    #donot check the size
    #(already checked at forward elimination)
    calc_result(c, [])
    |> modify_results_order(hist)
  end




  def solve_gaussian_eliminate(matrix, inhom_vector,
    pivoting \\ false, is_comp \\ false) do
    expanded_matrix = expand(matrix, inhom_vector)
    {ut_matrix, hist} = forward_eliminate(expanded_matrix, pivoting, is_comp)
    backward_substitute(ut_matrix, hist)
  end

end
