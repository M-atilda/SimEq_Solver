defmodule SimEq.Iteration do
  import Matrix

  # defolt params for iterative calcuration
  @defolt_tolerance 1.0e-06
  @defolt_max_iteration_times 1000

  #@fn     [float] solve_jacob(%Matrix, [float], ())
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
    product = Enum.zip(l, result)
              |> Enum.reduce(0, fn ({f, s}, acm) -> (f*s)+acm end)
    sub_val = product - (Enum.at(l, i-1) * Enum.at(result, i-1))
    (inhom_val - sub_val) / Enum.at(l, i-1)
  end
  defp calc_next_result(%Matrix{line: l_m} = matrix, inhom_vector,
                        result) do
    for i <- 1..l_m do i end
    |> Enum.map(
      &(Task.await(Task.async(fn -> calc_next_val(&1,
                    get_line(matrix, &1), Enum.at(inhom_vector, &1-1),
                    result) end))))
  end

  defp is_acceptable_result(tol, result, new_result) do
    Enum.zip(result, new_result)
    |> Enum.map(fn {f, s}
      -> Task.await(Task.async(fn -> (tol*f) > abs(f-s) end)) end)
    |> Enum.reduce(true, fn (b, acm) -> acm && b end)
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


  defp calc_next_result_sor(result, sample_new_result, w) do
    Enum.zip(result, sample_new_result)
    |> Enum.map(fn {f, s} -> f + (w * (s - f)) end)
  end

  defp sor(_, _,
           result,
           _,
           _,
           ite, max_ite) when ite > max_ite, do: {:over, ite, result}
  defp sor(matrix, inhom_vector,
           result,
           w, # relaxation factor
           tol,
           ite, max_ite) do
    sample_new_result = calc_next_result_gs(matrix, inhom_vector, result)
    new_result = calc_next_result_sor(result, sample_new_result, w)
    if is_acceptable_result(tol, result, new_result) do
      {:fini, ite, new_result}
    else
      sor(matrix, inhom_vector, new_result, w, tol, ite+1, max_ite)
    end
  end




def solve_jacob(matrix, inhom_vector,
                  start_x,
                  max_ite \\ @defolt_max_iteration_times,
                  tol \\ @defolt_tolerance) do
    {status, i, result} = jacob(matrix, inhom_vector, start_x, tol, 1, max_ite)
    if status == :fini do
      {i, result}
    else
      raise AlgorithmError,
        message: "Invalid algorithm error <solve_jacob@iteration.ex>
  tolerance: #{tol}
  iteration: #{max_ite}",
        value: result
      end
  end

  def solve_gauss_seidel(matrix, inhom_vector,
                         start_x,
                         max_ite \\ @defolt_max_iteration_times,
                         tol \\ @defolt_tolerance) do
    {status, i, result} = gauss_seidel(matrix, inhom_vector,
                                       start_x,
                                       tol,
                                       1, max_ite)
    if status == :fini do
      {i, result}
    else
      raise AlgorithmError,
        message: "Invalid algorithm error <solve_gauss_seidel@iteration.ex>
  tolerance: #{tol}
  iteration: #{max_ite}",
        value: result
      end
  end

  def solve_SOR(matrix, inhom_vector,
                start_x,
                relaxation_factor,
                max_ite \\ @defolt_max_iteration_times,
                tol \\ @defolt_tolerance) do
    {status, i, result} = sor(matrix, inhom_vector,
                              start_x,
                              relaxation_factor,
                              tol,
                              1, max_ite)
    if status == :fini do
      {i, result}
    else
      raise AlgorithmError,
        message: "Invalid algorithm error <solve_sor@iteration.ex>
  relaxation Factor: #{relaxation_factor}
  tolerance: #{tol}
  iteration: #{max_ite}",
        value: result
    end
  end


  #NOTE: this way is not the recommended way to get suitable relaxation factor
  #      it returns the value which satisfy the defolt tolerance fastest
  def test_relaxation_factor(matrix, inhom_vector,
                             start_x,
                             min, max, delta,
                             ite \\ @defolt_max_iteration_times) do
    next_w = Stream.iterate(min, &(&1+delta))
             |> Enum.take(round((max-min)/delta)+1)
             |> Enum.map(fn w
      -> {w, try do
               solve_SOR(matrix, inhom_vector,
                         start_x,
                         w,
                         ite,
                         @defolt_tolerance)
               true
             rescue
               _ ->
                 false
             end} end)
             |> Enum.reduce([], fn ({w, b}, acm)
               -> if b do
               [w|acm]
             else
               acm
             end end)
    if length(next_w) == 0 do
      if ite == @defolt_max_iteration_times do
        raise AlgorithmError,
          message: "in designated relaxation factor and defolt tolerance the equations cannot be solved"
      else
        #NOTE: this way is overbearing
        (min + max) / 2
      end
    else
      if length(next_w) == 1 do
        Enum.at(next_w, 0)
      else
          test_relaxation_factor(matrix, inhom_vector,
                                 start_x,
                                 Enum.min(next_w), Enum.max(next_w), delta,
                                 ite-1)
      end
    end
  end


  def compare_iteration(matrix, inhom_vector,
                        start_x,
                        w, # relaxation factor
                        max_ite,
                        tol) when is_number(max_ite) do
    {_, r_j} = try do
                 solve_jacob(matrix, inhom_vector,
                             start_x,
                             max_ite,
                             tol)
               rescue
                 ex ->
                   {0, AlgorithmError.get_value(ex)}
               end
    {_, r_gs} = try do
                  solve_gauss_seidel(matrix, inhom_vector,
                                     start_x,
                                     max_ite,
                                     tol)
                rescue
                  ex ->
                    {0, AlgorithmError.get_value(ex)}
                     end
    {_, r_sor} = try do
                   solve_SOR(matrix, inhom_vector,
                             start_x,
                             w,
                             max_ite,
                             tol)
                 rescue
                   ex ->
                     {0, AlgorithmError.get_value(ex)}
                 end
    [r_j, r_gs, r_sor]
    |> Enum.map(fn val -> SimEq.get_error_norm(matrix, inhom_vector, val) end)
  end
  def compare_residual(matrix, inhom_vector,
                       start_x,
                       relaxation_factor,
                       max_ite,
                       tol) do
    Enum.map(max_ite, fn ite -> compare_iteration(matrix, inhom_vector,
                                                  start_x,
                                                  relaxation_factor,
                                                  ite,
                                                  tol) end)
  end

end
