defmodule SimEq.Interface do

  def main(_argv) do
    #parse_args argv
    # coefficient = Matrix.generate([[5, -3, -1, 0, 2, 1],
    #                                [0, -5, 2, 1, 0, 2],
    #                                [-1, 3, -6, -2, 3, 3],
    #                                [-1, 0, 3, 5, -2, -1],
    #                                [0, 3, 3, -1, 5, -4],
    #                                [2, 3, 2, 3, 2, -4]])
    # inhomogenious = [16, 12, 12, 16, 18, 20]
    # start_x = [1, 1, 1, 1, 1, 1]


    # coefficient = Matrix.generate(
    #   [[2,0,0,0,-1,0,-0.5,0.5,0],
    #    [0,2,0,-1,0,0,0.5,-0.5,0],
    #    [0,0,2.5,0.5,-0.5,-0.5,-1,0,0],
    #    [0,-1,0.5,1.5,-0.5,-0.5,0,0,0],
    #    [-1,0,-0.5,-0.5,2,1,0,0,-0.5],
    #    [0,0,-0.5,-0.5,1,2,0,-1,0.5],
    #    [-0.5,0.5,-1,0,0,0,2.5,-0.5,-1],
    #    [0.5,-0.5,0,0,0,-1,-0.5,1.5,0],
    #    [0,0,0,0,-0.5,0.5,-1,0,1.5]])
    # p_3y_ex = -500 / (2.8 * :math.pow(10, 7))
    # p_5y_ex = -300 / (2.8 * :math.pow(10, 7))
    # # p_3y_ex = -500
    # # p_5y_ex = -300
    # inhomogenious = [0,0,0,p_3y_ex,0,0,0,p_5y_ex,0]
    # #start_x = [0,0,0,0,0,0,0,0,0]



    # IO.puts "@gaussian"
    # IO.puts "* result *"
    # result = SimEq.Gauss.solve_gaussian_eliminate(coefficient, inhomogenious)
    # Enum.map(result, fn val -> IO.puts val end)

    # IO.puts "* error *"
    # error = SimEq.get_error(coefficient, inhomogenious, result)
    # Enum.map(error, fn val -> IO.puts val end)
    # IO.puts "* error norm *"
    # IO.puts SimEq.get_error_norm(coefficient, inhomogenious, result)




    # IO.puts "@gauss-seidel"
    # IO.puts "* result *"
    # {ite, result} = SimEq.Iteration.solve_gauss_seidel(coefficient, inhomogenious, start_x, 1000, 1.0e-6)
    # IO.puts "iteration: #{ite}"
    # Enum.map(result, fn val -> IO.puts val end)

    # IO.puts "* error *"
    # error = SimEq.get_error(coefficient, inhomogenious, result)
    # Enum.map(error, fn val -> IO.puts val end)
    # IO.puts "* error norm *"
    # IO.puts SimEq.get_error_norm(coefficient, inhomogenious, result)



  # def parse_args(argv) do
  #   parse = OptionParser.parse(argv,
  #     swiches: [help: :boolean],
  #     aliases: [h: :help])

  #   # case parse do

  #   # end
  end

end
