# README

The parameters that may vary during testing (such as the mesh file path) have been extracted into an external `test_config.json` file (referred to as the `config file` below). This allows you to modify the test parameters without recompiling the test program.

## What is the config file?

The `config file` is organized in the format `TestSuitName.TestName.ParamName: Value`, corresponding to the two-level structure of `TestSuitName` and `TestName` in `Google Test`.
For each `Test` (uniquely identified by `TestSuitName` and `TestName`), you can define multiple parameters, such as:

```json
  "Boolean": {
    "TestIfCrash": {
      "mesh_path_1": "./data/boolean_data/bunny25k.obj",
      "mesh_path_2": "./data/boolean_data/cow.obj"
    }
  }
```

> * Of course, the `TestSuitName` and `TestName` in the `config file` do not need to exactly match those in `Google Test`.
> * You can modify the names in the `config file` as you see fit for your usage.
> * However, it is recommended that the names in the `config file` correspond to or resemble those in `Google Test` to facilitate searching and editing by others.

## How to use the config file?

* When running the test program, use the `--config=path-to-json` argument to allow the test program to read the `config file` (similar to `--gtest`...).
* To change the test parameters, simply modify the relevant values in the `config file`.
* File paths may differ on different test machines, so it is recommended to copy the `config file`, modify its contents as needed, and then adjust the test program's `--config` parameter to point to the new `config file`.

## How to make your Test utilize this feature

* First, add your `TestSuitName.TestName` structure in the `config file` and include your parameter structure `{Param:Value, ...}` (you can refer to existing settings). For example:

  ```json
  "FastQEM": {
    "Common": {
      "dir": "./data/",
      "filename": "bunny_34k.obj",
      "ratios": [
        0.01,
        0.1,
        0.3,
        0.6
      ]
    }
  }
  ```

* Next, in your `Test` function, retrieve these parameters. The macro `TEST_GET_CONFIG` is defined in the `test_utils.h` file and allows you to obtain the `config object` corresponding to a `Test`.
  * The macro takes the parameters `(TestSuitName, TestName)`, corresponding to the TestSuitName.TestName in the config file.
  * Call this macro at the beginning of your `Test` function to obtain the `config object` corresponding to `TestSuitName.TestName`.
  * For example:

    ```cpp
    TEST_GET_CONFIG(FastQEM, Common);
    ```

* The `config object` can be used as follows:

  * For parameters of types like `int`, `double`, `std::string`, etc., you can obtain the value using `config.get<ParamType>("ParamName")`. For instance, to get the value of `dir`:

    ```cpp
    std::string dir = config.get<std::string>("dir");
    ```

  * For array-type parameters, you can iterate through the array contents using an iterator. For example, you can use the following code in the Test function to get the values of `ratios`:

    ```cpp
    std::vector<double> ratios;
    for (auto it : config.get_child("ratios"))
      ratios.push_back(it.second.get_value<double>());
    ```
