defmodule AlgorithmError do
  defexception message: "Invalid algorithm is used",
    value: nil,
    can_retry: false

  def get_value(me) do
    me.value
  end
end
