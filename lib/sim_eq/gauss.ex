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




  def get_pivot_index(matrix, i, false) do
    #partial pivoting
    candidates = Enum.drop(get_row(matrix, i), i-1)
    pivot_num = max(Enum.max(candidates), abs(Enum.min(candidates)))
    {i + Enum.find_index(candidates, &(abs(&1) == pivot_num)), i}
  end
  def get_pivot_index(%Matrix{line: l_m} = matrix, i, true) do
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



  # helpers for forward elimination
  defp calc_divided_line({i, j, div_num}, matrix) do
    Task.await(Task.async(fn ->
          Enum.zip(get_line(matrix, i), get_line(matrix, j))
          |> Enum.map(fn {f, s} -> s-(f*div_num) end)
          # get_line(mul_R_l(matrix, i, j, -div_num), j)
        end))
  end

  #@fn     {%Matrix, [{int,int}]} fe_helper(%Matrix, int, boolean, [{int,int}])
  #@param  coefficient matrix while conversion to upper triangular one
  #@param  the number of lines already precessed
  #@param  use completion pivoting or not
  #@param  row swapping history
  defp fe_helper(%Matrix{line: l_m} = matrix,
    i, _, acm) when l_m == i, do: {matrix, acm}
  defp fe_helper(%Matrix{line: l_m, contents: c} = matrix, i,
    is_comp, acm) do
    {pivot_i, pivot_j} = get_pivot_index(matrix, i, is_comp)
    swaped_matrix = mul_P_l(matrix, i, pivot_i)
    |> mul_P_r(i, pivot_j)
    new_hist = if is_comp do
      [{i, pivot_j}|acm] else acm end
    pivot_num = get_factor(swaped_matrix, i, i)
    div_nums = for j <- (i+1)..l_m do
      {i, j, get_factor(swaped_matrix, j, i) / pivot_num} end
    new_matrix = generate(
      Enum.take(c, i) ++
      for x <- div_nums do
        calc_divided_line(x, swaped_matrix) end)
    fe_helper(new_matrix, i+1, is_comp, new_hist)
  end


  def forward_eliminate(%Matrix{line: l_m, row: r_n}, _)
    when l_m != r_n, do: raise "math error <forward_eliminate@gauss.ex>"
  def forward_eliminate(m, is_comp) do
    fe_helper(m, 1, is_comp, [])
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

  def modify_results_order(result, []), do: result
  def modify_results_order(result, [{f_id, s_id}|tail]) do
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
    is_comp \\ false) do
    expanded_matrix = expand(matrix, inhom_vector)
    {ut_matrix, hist} = forward_eliminate(expanded_matrix, is_comp)
    backward_substitute(ut_matrix, hist)
  end

end
