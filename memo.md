# 列名の変更
all_data.rename(columns={'text': '口コミ内容', 'label': '評価'}, inplace=True)

# インデックスを列に格納
# all_data = all_data.reset_index()

# フレーム作成
frame = all_data

# 左揃え設定
frame.style.set_properties(**{'text-align': 'left'})

# ファイル数分のループ処理
for i in range(0, 11):
    # 各ファイルの読み込み
    df = pd.read_csv(f'review_{i}.csv')
    # print(df.head())
    
# 列名変更
    df = df.rename(columns={'text': 'review_text', 'label': 'sentiment'})

    # データの結合
    if i == 1:
        all_data = df
    else:
        all_data = pd.concat([all_data, df], ignore_index=True)
    print(all_data)


# 前処理済みテキストの表示
print("前処理済みテキスト:")
print(all_data['review_text'])

# 前処理済みデータの保存
all_data.to_csv('preprocessed_data.csv', index=False)

# 不要な文字列の削除
all_data['review_text'] = all_data['review_text'].str.replace(r'[\[\]]', '')  # 角括弧
all_data['review_text'] = all_data['review_text'].str.replace('[\n]', '')  # 改行文字
all_data['review_text'] = all_data['review_text'].str.replace(r'[^\x00-\x7F]', '')  # 全角文字

# 特殊文字の処理
all_data['review_text'] = all_data['review_text'].str.replace(r'[。？！，；、・:：…―「」\'\']', '')  # 句読点

# 正規化
all_data['review_text'] = all_data['review_text'].str.lower()  # 小文字変換
all_data['review_text'] = all_data['review_text'].str.replace(' ', '_')  # 空白をアンダースコアに変換
all_data['review_text'] = all_data['review_text'].str.replace('　', '_')  # 全角空白をアンダースコアに変換
all_data['review_text'] = all_data['review_text'].str.replace('n', 'ん')  # ンをんに変換
all_data['review_text'] = all_data['review_text'].str.replace('っ', 'つ')  # ッをっに変換