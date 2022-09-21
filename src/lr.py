from bio import *
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import ElasticNet, LinearRegression

from src.consts import muted_columns
from src.helpers.model import prepare_model_data, get_accuracy, remove_outliers
from src.utils.common import inlogit, logit


def custom_score_grouped(groups):
    def custom_score(model, X, y):
        df = X.copy()
        df['fraq'] = inlogit(y)
        scores = get_accuracy(model, df, is_numpy=False)
        return 0 if np.isnan(scores['cor']) else scores['cor']
        # grouped = []
        # for _, group in groups:
        #     df = loc(X, group.index)
        #     df['fraq'] = inlogit(loc(y, group.index))
        #     scores = get_accuracy(model, df)
        #     grouped.append(0 if np.isnan(scores['uplift']) else 1 / (1 + abs(scores['uplift'])))
        #
        # return np.median(grouped)
    #
    return custom_score


def elastic_net(train, alpha=np.power(2.0, range(-5, 2)), l1_ratio=[1, 0.5, 0.1], **kwargs):
    train_ = remove_outliers(train)
    cols, train_X, train_Y = prepare_model_data(train_, is_numpy=False)
    if train_Y.empty or np.var(inlogit(train_Y)) < 0.001:
        model = LinearRegression()
        model.coef_ = np.zeros(shape=(len(cols), ))
        model.intercept_ = logit(np.median(train_['fraq']) if not train_.empty else train['fraq'].min())
        coefs = pd.DataFrame({'(Intercept)': {'Estimate': model.intercept_, 'p-value': 0}}).T
        return {
            'coefs': coefs,
            'cv': None,
            'model': model,
            'params': 'low variance',
        }

    custom_scorer = custom_score_grouped(train_.groupby(by='Tissue'))

    cv = GridSearchCV(
        estimator=ElasticNet(random_state=0),
        cv=2,
        param_grid={'alpha': alpha, 'l1_ratio': l1_ratio},
        scoring=custom_scorer,
    )
    cv.fit(train_X, train_Y, sample_weight=train_['Freq'] if 'Freq' in train_.columns else None)
    model = cv.best_estimator_
    significant_inds = model.coef_ != 0
    cols, train_X, train_Y = prepare_model_data(train_[set(list(np.array(cols)[significant_inds])+muted_columns)&set(train_.columns)], is_numpy=False)

    if not cols:
        coefs = pd.DataFrame({'(Intercept)': {'Estimate': model.intercept_, 'p-value': 0}}).T
        return {
            'coefs': coefs,
            'cv': None,
            'model': model,
            'params': 'low variance',
        }

    model = LinearRegression(**kwargs)
    model.fit(train_X, train_Y, sample_weight=train_['Freq'])

    coefs = pd.DataFrame(model.coef_, index=cols)
    coefs.loc[:, 'p-value'] = 0
    coefs.loc['(Intercept)'] = [model.intercept_, 0]
    coefs.columns = ['Estimate', 'p-value']

    return {
        'coefs': coefs,
        'cv': cv,
        'model': model,
        'params': cv.best_params_,
    }
