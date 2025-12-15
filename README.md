# Script Library  
&gt; 我的日常脚本库 —— 随写、随用、随优化

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Bash](https://img.shields.io/badge/lang-bash-%23239120)](https://www.gnu.org/software/bash/)
[![Python](https://img.shields.io/badge/lang-python-3776AB)](https://python.org)
[![GPT-powered](https://img.shields.io/badge/GPT-powered-74aa9c)](https://openai.com)

---

##  项目定位
本仓库收录了我在科研/服务器运维中**高频使用**的所有脚本，  
所有代码均通过 **ChatGPT 复盘+重构** 后同步至此，保证“即用即新”。  
&gt; 目标：**一次写好，终身复用，跨机零配置。**

---

##  脚本速查
| 目录 | 脚本 | 功能 | 更新日期 |
|------|------|------|----------|
| ./ | fasta_tool.py | 该脚本是一个功能全面的FASTA序列处理工具，支持按ID提取/剔除、长度过滤、去重、格式标准化、序列变换（反向互补、大小写、子序列）、碱基统计及合法性验证等 | 2025-12-14 |
| ./aa_cdr_gene_identification | run_cdr_genotype.py | 一键式抗体 CDR  Germline 分析工具（AbNumber 提取 CDR, ANARCI 基因型分配,结果合并） | 2025-12-13 |
| ./ | dusub_submission.sh | 批量通过 nohup 提交 '<scripts_dir>' 目录下的所有 '.sh' 脚本 | 2025-12-15 |
| ./server_admin_script | create_bio_user_lock_passwd.sh | 创建生信用户并禁止其自行修改密码 | 2025-11-15 |

 **完整地图**：打开 [`SCRIPT_MAP.md`](SCRIPT_MAP.md) 获取参数说明、依赖与示例。

---

##  快速开始
```bash
# 1. 克隆即得
git clone https://github.com/AdolfFK/script_library.git
cd script_library

# 2. 把 bin 加入全局 PATH（可选）
echo "export PATH=\"$PWD/bin:\$PATH\"" &gt;&gt; ~/.bashrc
source ~/.bashrc

# 3. 直接开跑
dusub_submission.sh ./my_jobs