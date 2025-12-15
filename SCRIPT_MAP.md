| 目录 | 脚本 | 功能 | 更新日期 |
|------|------|------|----------|
| ./ | fasta_tool.py | 该脚本是一个功能全面的FASTA序列处理工具，支持按ID提取/剔除、长度过滤、去重、格式标准化、序列变换（反向互补、大小写、子序列）、碱基统计及合法性验证等 | 2025-12-14 |
| ./aa_cdr_gene_identification | run_cdr_genotype.py | 一键式抗体 CDR  Germline 分析工具（AbNumber 提取 CDR, ANARCI 基因型分配,结果合并） | 2025-12-13 |
| ./ | dusub_submission.sh | 批量通过 nohup 提交 '<scripts_dir>' 目录下的所有 '.sh' 脚本 | 2025-12-15 |
| ./server_admin_script | create_bio_user_lock_passwd.sh | 创建生信用户并禁止其自行修改密码 | 2025-11-15 |
| ./server_admin_script | create_bio_user_change_passwd.sh | 创建生信用户并初始化工作目录 | 2025-12-15 |
| ./Single_cell_sequence | run_MultiSamples.V5_gpt_v3.R | 单细胞的第一步脚本，处理cellranger的结果文件，可以做分组。结果最好还是人工审核一下。把图上的点变的更小 | 2025-12-15 |
| ./Single_cell_sequence | generate_batch_script.py | 这个脚本用来批量生成单细胞初步的脚本 | 2025-06-10 |
| ./Single_cell_sequence | run_MultiSamples.V5_gpt_v2.R | 单细胞的第一步脚本，处理cellranger的结果文件，可以做分组。结果最好还是人工审核一下 | 2025-12-15 |
| ./ | map_generator.sh | 递归扫描子目录，生成 SCRIPT_MAP.md（含目录列） | 2025-12-14 |
| ./ | sort_fasta_universal.py | 对 FASTA 文件多维度排序（ID 字母、ID 数字、序列长度、正/反序） | 2025-12-10 |
| ./ | sort_fasta_by_id.py | 按序列ID对FASTA文件进行排序，并按标准格式输出,默认每行60字符 | 2025-12-09 |
| ./ | sync2gh.sh | 极简版——永远向 main 分支提交并推送；版本管理完全由用户自行处理。 | 2025-12-13 |
| ./trans2cmp | do_compare.R | compare_enrichment.R脚本的运行样例 | 2025-12-14 |
| ./trans2cmp | compare_enrichment.R | 脚本比较两次转录组分析的结果文件，输出其中的差异 | 2025-12-14 |
