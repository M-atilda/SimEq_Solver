defmodule SimEq.Iteration do
  import SimEq.Matrix

  # defolt params for iterative calcuration
  @defolt_tolerance 1.0e-6
  @defolt_max_iteration_times 100

  #@fn     [float] solve_jacob(%SimEq.Matrix, [float], ())
  #@brief  this function returns the taple
  #          which consists of the end-status,
  #                            the times of iteration
  #                            and result vector
  #        end-status is represented as a simbol like below
  #          :fini   the error size becomes smaller than designated tolerance
  #          :over   reach the max iteration times

  @moduledoc """
  provides iterative methods
    for solving simultenious equations

  the calcuration finishes
    when all the value of the result are
      within the margen of error
  the tolerance is measured at each factor
  """

  #TODO:
  #  * there are some duplicated codes
  #    make generic function like facade pattern




  defp calc_next_val(i,
                     l, inhom_val,
                     result) do
    sub_val = Enum.zip(l, result)
                |> List.delete({Enum.at(l, i-1), Enum.at(result, i-1)})
                |> Enum.reduce(0, fn ({f, s}, acm) -> (f*s)+acm end)
    (inhom_val - sub_val) /
    Enum.at(l, i-1)
  end
  defp calc_next_result(%SimEq.Matrix{line: l_m} = matrix, inhom_vector,
                        result) do
    for i <- 1..l_m do i end
    |> Enum.map(fn i
      -> calc_next_val(i,
                       get_line(matrix, i), Enum.at(inhom_vector, i-1),
                       result) end)
  end

  defp is_acceptable_result(tol, result, new_result) do
    Enum.zip(result, new_result)
    |> Enum.reduce(true, fn ({f, s}, acm) -> acm && (tol > abs(f-s)) end)
  end

  defp jacob(_, _,
             result,
             _,
             ite, max_ite)
  when ite > max_ite, do: {:over, ite, result}
  defp jacob(matrix, inhom_vector,
             result,
             tol,
             ite, max_ite) do
    new_result = calc_next_result(matrix, inhom_vector, result)
    if is_acceptable_result(tol, result, new_result) do
      {:fini, ite, new_result}
    else
      jacob(matrix, inhom_vector, new_result, tol, ite+1, max_ite)
    end
  end


  defp cnr_gs_helper(_,
                     _, _,
                     result, []), do: result
  defp cnr_gs_helper(i,
                     matrix, [inhom_val|inhom_tail],
                     new_result, result) do
    next_val = calc_next_val(i, get_line(matrix, i),
      inhom_val, new_result++result)
    cnr_gs_helper(i+1,
      matrix, inhom_tail,
      new_result++[next_val], Enum.drop(result, 1))
  end
  defp calc_next_result_gs(matrix, inhom_vector,
                           result) do
    cnr_gs_helper(1, matrix, inhom_vector, [], result)
  end

  defp gauss_seidel(_, _,
                    result,
                    _,
                    ite, max_ite) when ite > max_ite, do: {:over, ite, result}
  defp gauss_seidel(matrix, inhom_vector,
                    result,
                    tol,
                    ite, max_ite) do
    new_result = calc_next_result_gs(matrix, inhom_vector, result)
    if is_acceptable_result(tol, result, new_result) do
      {:fini, ite, new_result}
    else
      gauss_seidel(matrix, inhom_vector, new_result, tol, ite+1, max_ite)
    end

  end




  def solve_jacob(matrix, inhom_vector,
                  start_x,
                  tol \\ @defolt_tolerance,
                  max_ite \\ @defolt_max_iteration_times) do
    {status, i, result} = jacob(matrix, inhom_vector, start_x, tol, 1, max_ite)
    if status == :fini do
      {i, result}
    else
      raise AlgorithmError,
          message: "Invalid algorithm error <solve_jacob@iteration.ex>
  start value: #{start_x}
  tolerance: #{tol}
  iteration: #{max_ite}"
    end
  end

  def solve_gauss_seidel(matrix, inhom_vector,
                         start_x,
                         tol \\ @defolt_tolerance,
                         max_ite \\ @defolt_max_iteration_times) do
    {status, i, result} = gauss_seidel(matrix, inhom_vector,
      start_x,
      tol,
      1, max_ite)
    if status == :fini do
      {i, result}
    else
      raise AlgorithmError,
          message: "Invalid algorithm error <solve_gauss_seidel@iteration.ex>
  start value: #{start_x}
  tolerance: #{tol}
  iteration: #{max_ite}"
    end
  end

end
