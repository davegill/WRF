# Jenkins Testing Directory

For users looking to modify the testing of WRF, look in the `terraform` directory, the `main.tf` file. This has the list of all tests, explicitly listed. A user may reduce the number of tests by modifying the value of `default = 10` to something like `default = 1`.
```
//Run number of test case based on count value: say 3 will run 3 test cases 
variable "instance_count" {
 default = 10
}
```

