defmodule SimEq.Gauss do
  import SimEq.Matrix

  #@fn      [float] solve_gaussian_eliminate(%SimEq.Matrix, [float], (bool))
  #@params  coefficient matrix
  #@params  inhomogenious vector
  #@params  use pivoting or not(option)

  #NOTE:
  #  * the history of swapping lines
  #        for use suitable value as pivot
  #      is memorized as taples' list
  #    ex) swap line 1 and 3 => [{1, 3}]
  #
  #  * pivot is choiced in partial pivoting
  #    complete pivoting is not yet implemeted
  #
  #  * whether the simultinious equations is solvable
  #      is not treated at this module
  #    it assumes pre-process has already been done
  #      (its rank is same as the number of variables)




  # helpers for forward elimination
  defp get_pivot_index(matrix, i) do
    #partial pivoting
    targets_list = Enum.drop(get_row(matrix, i), i-1)
    pivot_num = max(Enum.max(targets_list), abs(Enum.min(targets_list)))
    i + Enum.find_index(targets_list, &(abs(&1) == pivot_num))
  end

  defp div_sub_lines([], _, _, matrix) do: matrix
  defp div_sub_lines([h|tail], j, i, matrix) do
    div_sub_lines(tail, j+1, i, mul_R_l(matrix, j, i, -h))
  end

  defp fe_inner(%SimEq.Matrix{line: l_m} = matrix, i,
    _, acm \\ []) when l_m == i do: {matrix, acm}
  defp fe_inner(%SimEq.Matrix{line: l_m, contents: c} = matrix, i,
    false, _ \\ nil) do
    pivot_num = get_factor(matrix, i, i)
    div_nums = for j <- (i+1)..l_m, do: get_factor(Matrix, j, i) / pivot_num
    fe_inner(div_sub_lines(div_nums, i+1, i, matrix),
      i+1,
      false)
  end
  defp fe_inner(%SimEq.Matrix{line: l_m, contents: c} = matrix, i,
    true, acm) do
    pivot_id = get_pivot_index(matrix, i)
    swaped_matrix = mul_P_l(matrix, i, pivot_id)
    pivot_num = get_factor(matrix, i, i)
    div_nums = for j <- (i+1)..l_m, do: get_factor(Matrix, j, i) / pivot_num
    fe_inner(div_sub_lines(div_nums, i+1, i, matrix),
      i+1,
      true,
      [{i, pivot_id}|acm])
  end


  def forward_eliminate(%SimEq.Matrix{line: l_m, row: r_n}, _)
    when l_m != r_n do: raise "math error <forward_eliminate@gauss.ex>"
  def forward_eliminate(m, pivoting) do: fe_inner(m, 1, pivoting)



  # helpers for backward substitution
  defp calc_result([], result) do: result
  defp calc_result(c, []) do
    calc_result(Enum.slice(c, 0..-2), [List.last(List.last c)])
  end
  defp calc_result(c, result) do
    [inhom_val, last_l] = Enum.reverse(List.last c)
    sub_val = Enum.zip(Enum.reverse(result), last_l)
           |> Enum.reduce fn ({f, s}, acm) -> (f*s)+acm end
    calc_result(Enum.slice(c, 0..-2), [inhom_val-sub_val|result])
  end

  defp modify_results_order(result, []) do: result
  defp modify_results_order(result, [{f_id, s_id}|tail]) do
    f_val = Enum.at(result, f_id)
    s_val = Enum.at(result, s_id)
    new_result = List.replace_at(result, f_id-1, s_val)
              |> List.replace_at(s_id-1, f_val)
    modify_results_order(new_result, tail)
  end


 s def backward_substitute(%SimEq.Matrix{contents: c}, hist) do
    #donot check the size
    #(already checked at forward elimination)
    result = calc_result(c, [])
          |> modify_results_order hist
  end




  def solve_gaussian_eliminate(matrix, inhom_vector, pioting \\ true) do
    expanded_matrix = expand(matrix, inhom_vector)
    {ut_matrix, hist} = forward_elimination(expanded_matrix, pivoting)
    backward_substitute(ut_matrix, hist)
  end

end
