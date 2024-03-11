# README

测试中可变的参数（如网格文件路径）已经提取到外部的test_config.json文件（下面简称config文件）中，方便修改测试参数，且不需要重新编译测试程序。

## config文件是怎样的？

`config文件`的组织形式为`TestSuitName.TestName.ParamName: Value`，对应`google test`中的两级组织形式`TestSuitName`和`TestName`。  
对每个`Test`（由`TestSuitName+TestName`唯一确定），可以定义多个参数，如

```json
  "Boolean": {
    "TestIfCrash": {
      "mesh_path_1": "./data/boolean_data/bunny25k.obj",
      "mesh_path_2": "./data/boolean_data/cow.obj"
    }
  }
```

> 当然，`config文件`中的`TestSuitName`和`TestName`并没有限制要完全和`google test`中的对应一致。  
> 你可以任意修改`config文件`中的名字，以方便你的使用。  
> 建议`config文件`和`google test`中的名称对应或相似，以方便其他人查找和修改。

## 如何使用config文件？

* 在运行测试程序时，输入参数`--config=path-to-json`让测试程序可以读取`config文件`（和`--gtest...`类似）。
* 要改变测试参数时，改变`config文件`中的相关值即可。
* 不同测试机器上的文件路径可能不一样，请自行复制一份`config文件`修改内容，并修改测试程序的`--config`参数到新的`config文件`。

## 如何让你写的Test利用上这个特性

* 首先，在`config文件`中添加你的`TestSuitName.TestName`结构体，并添加你的参数结构体`{Param:Value, ...}`（可参考已存在的设置）。比如

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

* 其次，在你的`Test`函数中获取这些参数。`test_utils.h`文件中定义了一个宏`TEST_GET_CONFIG`，通过该宏得到`Test`对应的`config对象`
  * 该宏的参数是`(TestSuitName, TestName)`，对应`config文件`中的`TestSuitName.TestName`
  * 在你的`Test`函数开头调用该宏，来得到对应`TestSuitName.TestName`的`config对象`
  * 比如

    ```cpp
    TEST_GET_CONFIG(FastQEM, Common);
    ```

* `config对象`的使用方法是
  * 对于`int`、`double`、`std::string`等参数类型，可以通过`config.get<ParamType>("ParamName")`获得参数对应的值。比如获取`dir`的值

    ```cpp
    std::string dir = config.get<std::string>("dir");
    ```

  * 对于数组类型的参数，可以通过迭代器遍历数组内容。比如，可以在Test函数中使用如下代码获取`ratios`的值

    ```cpp
    std::vector<double> ratios;
    for (auto it : config.get_child("ratios"))
      ratios.push_back(it.second.get_value<double>());
    ```
